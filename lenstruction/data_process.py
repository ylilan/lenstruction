from __future__ import print_function
__author__ = 'lilan yang'

from lenstronomy.Util import kernel_util
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clipped_stats,SigmaClip
from astropy import wcs
from photutils import detect_threshold, detect_sources,deblend_sources, source_properties, Background2D, SExtractorBackground,make_source_mask
from photutils.datasets import make_noise_image
from six.moves import input
from lenstronomy.Data.imaging_data import ImageData
#from galight.data_process import DataProcess as galightDataProcess



class DataProcess(object):
    """
    The class contains useful fuctions to do, e.g. cut image, calculate cutsize, make mask of images.....
    """
    def __init__(self, image_data, header, snr=3.0, npixels=20, deltaPix=None, exp_time=None,background_rms=None,
                 wcs_data=None, background=None, segmap=None, psf=None, kernel = None, interaction = True):
        """

        :param data: fits format, image fits file.
        :param snr: float, signal-to-noise value.
        :param npixels: int, number of connected pixels that can be detected as a source.
        :param deltaPix: float, pixel size of the fits files.
        :param exp_time:  float, exposure time of the fits files.
        :param background_rms: float,float or array_like, the gaussian 1-sigma background noise in data.
        :param background:  float or 2D array,background value of the input image.
        :param kernel: The 2D array, filter the image before thresholding.
        :param interaction: bool, default is True.
        """

        self.header = header

        if wcs_data is None:
            self.wcs = wcs.WCS(header)
        else:
            print ("wcs is not none ")
            self.wcs = wcs_data
        if deltaPix is None:
            self.deltaPix= round((self.wcs.wcs_pix2world([[0, 100]], 1)[0][1]
                                  - self.wcs.wcs_pix2world([[0, 0]], 1)[0][1]) * 3600 * 0.01, 2)
        else:
            self.deltaPix = deltaPix
        self.snr = snr
        self.npixels = npixels
        if exp_time is None:
            print ("Note: Lenstruction read 'EXPTIME' keyword in the header !")
            exp_time = header['EXPTIME']
        else:
            exp_time = exp_time
        self.exp_time = exp_time
        self.background_rms_input = background_rms
        self.background_rms = None
        self.bakground = background
        self.kernel = kernel
        self.interaction= interaction
        self.image = image_data#self.hdul[0].data
        self.psf_input = psf
        if segmap is not None:
            self.segmap = fits.open(segmap)
        else:
            self.segmap = None

    def radec2detector(self,ra,dec):
        """
        WCS transformation from world to pixel coordinates
        :param ra: ra in deg
        :param dec: dec in deg
        :return:
        """
        y_detector = np.int(self.wcs.wcs_world2pix([[ra, dec]], 1)[0][0])
        x_detector = np.int(self.wcs.wcs_world2pix([[ra, dec]], 1)[0][1])

        return x_detector, y_detector

    def bkg_rms(self,x,y, r_cut):
        """
        :param background_rms:
        :return:
        """
        background_rms = self.background_rms_input

        if background_rms is not None:
            if type(background_rms) == type(1.0):
                background_rms_map =np.zeros((r_cut*2+1, r_cut*2+1)) + background_rms
            else:
                background_rms_map = fits.open(background_rms)[0].data[x - r_cut:x + r_cut + 1, y - r_cut:y + r_cut + 1]
        else:
            background_rms_map = background_rms

        self.background_rms = background_rms_map
        return background_rms_map

    def sub_bkg(self,img):
        sigma_clip = SigmaClip(sigma=3.)
        bkg_estimator = SExtractorBackground()
        mask_0 = make_source_mask(img, nsigma=3, npixels=5, dilate_size=11)
        mask_1 = (np.isnan(img))
        mask = mask_0 + mask_1
        bkg = Background2D(img, (5, 5), filter_size=(3, 3), sigma_clip=sigma_clip,
                           bkg_estimator=bkg_estimator, mask=mask)
        back = bkg.background * ~mask_1
        return img - back
        #return img

    def cut_image(self, x, y, r_cut):
        """
         Function used to cut input image.
         :param x: int, x coordinate in pixel unit
         :param y: int, y coordinate in pixel unit
         :param r_cut: int format value, radius of cut out image
         :return: cutted image
         """
        image_cutted_raw = self.image[x - r_cut:x + r_cut + 1, y - r_cut:y + r_cut + 1]
        image_cutted = self.sub_bkg(image_cutted_raw)
        return image_cutted

    def cut_image_psf(self, x, y, r_cut):
        """
         Function used to cut input image.
         :param x: int, x coordinate in pixel unit
         :param y: int, y coordinate in pixel unit
         :param r_cut: int format value, radius of cut out image
         :return: cutted image
         
         """
        image_cutted = self.image[x - r_cut:x + r_cut + 1, y - r_cut:y + r_cut + 1]
        return image_cutted

    def _seg_image(self, x, y, r_cut=100):
        """
        detect and deblend sources into segmentation maps
        :param x: int, x coordinate in pixel unit
        :param y: int, y coordinate in pixel unit
        :param r_cut:  int format value, radius of cut out image
        :return:
        """
        snr=self.snr
        npixels=self.npixels
        bakground = self.bakground
        error= self.bkg_rms(x,y,r_cut)
        kernel = self.kernel
        image_cutted = self.cut_image(x,y,r_cut)
        image_data = image_cutted
        threshold_detect_objs=detect_threshold(data=image_data, nsigma=snr,error=error)
        segments=detect_sources(image_data, threshold_detect_objs, npixels=npixels, filter_kernel=kernel)
        segments_deblend = deblend_sources(image_data, segments, npixels=npixels,nlevels=10)
        segments_deblend_info = source_properties(image_data, segments_deblend)
        nobjs = segments_deblend_info.to_table(columns=['id'])['id'].max()
        xcenter = segments_deblend_info.to_table(columns=['xcentroid'])['xcentroid'].value
        ycenter = segments_deblend_info.to_table(columns=['ycentroid'])['ycentroid'].value
        image_data_size = np.int((image_data.shape[0] + 1) / 2.)
        dist = ((xcenter - image_data_size) ** 2 + (ycenter - image_data_size) ** 2) ** 0.5
        c_index = np.where(dist == dist.min())[0][0]
        center_mask=(segments_deblend.data==c_index+1)*1 #supposed to be the data mask
        obj_masks = []
        for i in range(nobjs):
            mask = ((segments_deblend.data==i+1)*1)
            obj_masks.append(mask)
        xmin = segments_deblend_info.to_table(columns=['bbox_xmin'])['bbox_xmin'].value
        xmax = segments_deblend_info.to_table(columns=['bbox_xmax'])['bbox_xmax'].value
        ymin = segments_deblend_info.to_table(columns=['bbox_ymin'])['bbox_ymin'].value
        ymax = segments_deblend_info.to_table(columns=['bbox_ymax'])['bbox_ymax'].value
        xmin_c, xmax_c = xmin[c_index], xmax[c_index]
        ymin_c, ymax_c = ymin[c_index], ymax[c_index]
        xsize_c = xmax_c - xmin_c
        ysize_c = ymax_c - ymin_c
        if xsize_c > ysize_c:
            r_center = np.int(xsize_c)
        else:
            r_center = np.int(ysize_c)
        center_mask_info= [center_mask, r_center, xcenter, ycenter, c_index]
        return obj_masks, center_mask_info, segments_deblend


    def cutsize(self,x,y,r_cut=100,font_size=20):
     """
      determine cutsize of the image
     :param x: int, x coordinate in pixel unit
     :param y: y int, coordinate n pixel unit
     :param r_cut: int, radius of cut out image
     :return: cutout size
     """
     cutsize_data = r_cut
     if self.interaction & False:
         yn = True
         while(yn):
            m_image = self.cut_image(x, y, r_cut)
            fig_ci=plt.figure()
            plt.imshow(np.log10(m_image), origin='lower',cmap="gist_heat")
            plt.title('Good frame size? ('+repr(cutsize_data*2+1)+'x'+repr(cutsize_data*2+1)+' pixels^2' + ')',fontsize=font_size)
            plt.show(fig_ci)
            cutyn = input('Hint: appropriate frame size? (y/n): ')
            if cutyn == 'n':
                cutsize_ = np.int(input('Hint: input an appropriate cutsize (int format)?  (framesize=2*cutsize+1): '))
                cutsize_data = cutsize_
                r_cut = cutsize_
            elif cutyn == 'y':
                yn = False
                cutsize_data=cutsize_data
            else:
                raise ValueError("Please input 'y' or 'n' !")
     return cutsize_data



    def data_assemble(self, x,y, r_cut, add_mask=5,  pick_choice=False):
       """
       Function to pick up the pieces of source, lens light, pollution(mask)
       :param x: x coordinate in pixel unit
       :param y: y coordinate in pixel unit
       :param r_cut: radius size of the data.
       :param add_mask: number of pixels adding around picked pieces
       :param pick_choice: True to only select source region to fit, False to use whole observation to fit
       :return: kwargs_data
       """
       #segmentation components
       obj_masks,center_mask_info, segments_deblend_list = self._seg_image(x, y, r_cut=r_cut)
       data_masks_center, _, xcenter, ycenter, c_index = center_mask_info
       image = self.cut_image(x,y,r_cut)
       self.raw_image = image
       src_mask = np.zeros_like(image)
       lens_mask = np.zeros_like(image)
       plu_mask = np.zeros_like(image)
       lenslight_mask_index = []
       if self.segmap is not None and self.interaction:
           segmap=self.segmap[0].data
           segdata = segmap[x - r_cut:x + r_cut + 1, y - r_cut:y + r_cut + 1]
           plt.imshow(segdata, origin='lower')
           nlabel = np.unique(segdata)
           for i in range(nlabel.shape[0] - 1):
               ax = (int((np.where(segdata == nlabel[i + 1])[0].max() - np.where(segdata == nlabel[i + 1])[0].min()) / 2 +
                   np.where(segdata == nlabel[i + 1])[0].min()))
               ay = (int((np.where(segdata == nlabel[i + 1])[1].max() - np.where(segdata == nlabel[i + 1])[1].min()) / 3 +
                   np.where(segdata == nlabel[i + 1])[1].min()))
               plt.text(ay, ax, repr(nlabel[i + 1]), color='r', fontsize=15)
           plt.title('Input segmentation map')
           plt.show()
           source_mask_index = [int(sidex) for sidex in input('Selection of data via (inputed) segmentation index separated by space, e.g., 0 1 :').split()]
           for i in source_mask_index:
               src_mask = src_mask + segdata*(segdata==i*1)
        # lens light
           lenslightyn = input('Hint: is there lens light? (y/n): ')
           if lenslightyn == 'y':
                lenslight_mask_index = [int(lidex) for lidex in input('Selection of lens-plane light via (inputed) segmentation index separated by space, e.g., 0 1 :').split()]
                for i in lenslight_mask_index:
                    lens_mask = (lens_mask + segdata*(segdata==i*1))
           elif lenslightyn == 'n':
                lenslight_mask_index = []
           else:
                raise ValueError("Please input 'y' or 'n' !")
       # contamination
           pluyn = input('Hint: is there contamination? (y/n): ')
           if pluyn == 'y':
                plution_mask_index = [int(pidex) for pidex in input('Selection of contamination via (inputed) segmentation index separated by space, e.g., 0 1 :').split()]
                for i in plution_mask_index:
                    plu_mask = (plu_mask + segdata*(segdata==i*1))
           elif pluyn == 'n':
                plu_mask = np.zeros_like(image)
           else:
                raise ValueError("Please input 'y' or 'n' !")



       if self.segmap is None and self.interaction:
            self.plot_segmentation(image, segments_deblend_list, xcenter, ycenter, c_index)
            #source light
            if pick_choice:
                source_mask_index = [int(sidex) for sidex in input('Selection of data via segmentation index separated by space, e.g., 0 1 :').split()]
                for i in source_mask_index:
                    src_mask = src_mask + obj_masks[i]
            #lens light
            lenslightyn = input('Hint: is there lens light? (y/n): ')
            if lenslightyn == 'y':
                lenslight_mask_index = [int(lidex) for lidex in input('Selection of lens-plane light via segmentation index separated by space, e.g., 0 1 :').split()]
                for i in lenslight_mask_index:
                    lens_mask = (lens_mask + obj_masks[i])
            elif lenslightyn == 'n':
                lenslight_mask_index = []
            else:
                raise ValueError("Please input 'y' or 'n' !")
            # contamination
            pluyn = input('Hint: is there contamination? (y/n): ')
            if pluyn == 'y':
                plution_mask_index = [int(pidex) for pidex in input('Selection of contamination via segmentation index separated by space, e.g., 0 1 :').split()]
                for i in plution_mask_index:
                    plu_mask = (plu_mask + obj_masks[i])
            elif pluyn == 'n':
                    plu_mask =  np.zeros_like(image)
            else:
                raise ValueError("Please input 'y' or 'n' !")
       else:
            src_mask = data_masks_center


       #adding pixels around the selected masks
       selem = np.ones((add_mask, add_mask))
       src_mask = ndimage.binary_dilation(src_mask.astype(np.bool), selem)
       plu_mask_out = ndimage.binary_dilation(plu_mask.astype(np.bool), selem)
       plu_mask_out = (plu_mask_out - 1)*-1

       #select source region to fit, or to use whole observation to fit
       ##1.select source region to fit
       snr = self.snr
       source_mask = image * src_mask
       #create background image for picked
       if self.background_rms is None:
            _, _, std = sigma_clipped_stats(image, sigma=snr, mask=source_mask)
            tshape = image.shape
            img_bkg = make_noise_image(tshape, distribution='gaussian', mean=0., stddev=std, seed=12)
       else:
           tshape = image.shape
           std=np.mean(self.background_rms)
           img_bkg = make_noise_image(tshape, distribution='gaussian', mean=0., stddev=std, seed=12)

       no_source_mask = (src_mask * -1 + 1) * img_bkg
       picked_data = source_mask + no_source_mask

       ##2.use whole observation to fit while mask out the contamination
       maskedimg = image * plu_mask_out

       ##orginize the output 'kwargs_data'
       kwargs_data = {}
       if pick_choice:
           kwargs_data['image_data'] = picked_data#select source region to fit
       else:
           kwargs_data['image_data'] = maskedimg#use whole observation to fit while mask out the contamination

       if self.background_rms is None:
            kwargs_data['background_rms'] = std
            self.background_rms = std
       else:
            kwargs_data['background_rms'] = np.mean(self.background_rms)
       kwargs_data['exposure_time'] = self.exp_time
       kwargs_data['transform_pix2angle'] = np.array([[1, 0], [0, 1]]) * self.deltaPix
       ra_at_xy_0 = (y - r_cut) * self.deltaPix  # (ra,dec) is (y_img,x_img)
       dec_at_xy_0 = (x - r_cut) * self.deltaPix
       kwargs_data['ra_at_xy_0'] = ra_at_xy_0
       kwargs_data['dec_at_xy_0'] = dec_at_xy_0

       #coordinate of the lens light
       xlenlight, ylenlight = [], []
       if lenslight_mask_index !=[]:
            for i in lenslight_mask_index:
                xlenlight.append(ra_at_xy_0 + int(xcenter[i]) * self.deltaPix )
                ylenlight.append(dec_at_xy_0 + int(ycenter[i])* self.deltaPix )

       #for output
       self.data = kwargs_data['image_data']
       self.kwargs_data = kwargs_data
       self.data_mask = src_mask
       self.lens_mask  = lens_mask
       self.plu_mask = plu_mask_out
       self.obj_masks = obj_masks
       imageData = ImageData(**kwargs_data)
       self.imageData = imageData
       kwargs_seg = [segments_deblend_list, xcenter, ycenter, c_index]

       return kwargs_data, kwargs_seg, [xlenlight, ylenlight]


    def pick_psf(self, ra=None, dec=None, r_cut=50, pixel_size=None, kernel_size=None, pick=True):
        """
        select psf, either by coordinate or from input directly
        :param x:  x coordinate.
        :param y:  y coordinate.
        :param r_cut: radius size of the psf.
        :param deltaPix: pixel size of the psf.
        :param kernel_size: kernel size of the psf.
        :return: kwargs_psf
        """
        if self.psf_input is not None:
            image_psf = self.psf_input
        elif ra is not None:
            x_psf, y_psf = self.radec2detector(ra, dec)
            image_psf = self.cut_image_psf(x_psf, y_psf, r_cut)
        else:
            #automatically select psf in the field
            #data_process = galightDataProcess(fov_image=self.image, target_pos = [0,0], pos_type = 'pixel',                                                   )
            # print ("try to find psf by code itself")
            # if self.interaction & pick:
            #     psfyn = input('Hint: do you want to pick up the PSF youself? (y/n): ')
            #     if psfyn =='y':
            #         print(('Please only pick 1 psf, otherwise only the 1st will be chosen'))
            #         data_process.find_PSF(radius=r_cut, user_option=True)
            #         image_psf = data_process.PSF_list[0]
            #     elif psfyn =='n':
            #         data_process.find_PSF(radius=r_cut, user_option=False)
            #         image_psf = data_process.PSF_list[0]
            #     else:
            #         raise ValueError("Please input 'y' or 'n' !")
            # else:
            #     data_process.find_PSF(radius=r_cut, user_option=False)
            #     image_psf = data_process.PSF_list[0]
            print ("Try add psf function from Galight")
            # TODO Try add psf function from Galight
        if kernel_size is None:
            kernel_size = np.shape(image_psf)[0]

        image_psf_cut = kernel_util.cut_psf(image_psf, psf_size=kernel_size)#re-sized and re-normalized PSF


        if pixel_size is None:
            pixel_size=self.deltaPix
        else:
            pixel_size=pixel_size

        kwargs_psf = {'psf_type': 'PIXEL', 'kernel_point_source': image_psf_cut, 'pixel_size': pixel_size}
        self.psf = image_psf_cut
        return kwargs_psf



    def params(self, ra, dec, ra_psf = None, dec_psf =None, r_cut=100, add_mask=5, pick_choice=False,
               multi_band_type='joint-linear', kwargs_numerics={}, img_name='prodata_psf.pdf',img_id=0, if_plot=True):
        """
         image data parameters configuration in lenstronomy keywords arguments
        :param ra:
        :param dec:
        :param ra_psf:
        :param dec_psf:
        :param r_cut:
        :param add_mask:
        :param multi_band_type:
        :return: x, y coordinate in pixel units, image data keywords arguments
        """
        kwargs_numerics = self.numerics(**kwargs_numerics)
        kwargs_psf = self.pick_psf(ra=ra_psf, dec=dec_psf)
        multi_band_list =[]
        x_detector = []
        y_detector = []
        xylenslight_list = []
        data_mask_list = []
        lens_mask_list = []
        plu_mask_out_list = []
        for i in range(len(ra)):
            xy = self.radec2detector(ra[i], dec[i])
            x_detector.append(xy[0])
            y_detector.append(xy[1])
            cutsize = self.cutsize(xy[0], xy[1], r_cut=r_cut)
            kwargs_data, _, xylenslight= self.data_assemble(x=xy[0], y=xy[1],r_cut=cutsize, add_mask=add_mask,pick_choice=pick_choice)
            multi_band_list.append([kwargs_data, kwargs_psf, kwargs_numerics])
            xylenslight_list.append(xylenslight)
            data_mask_list.append(self.data_mask)
            lens_mask_list.append(self.lens_mask)
            plu_mask_out_list.append(self.plu_mask)
            #plot data, lens light, pollution
            if if_plot:
                self.plot_prodata_psf(img_name=img_name,img_id=img_id)
        kwargs_data_joint = {'multi_band_list': multi_band_list, 'multi_band_type': multi_band_type}
        self.data_mask_list = data_mask_list
        self.lens_mask_list = lens_mask_list
        self.plu_mask_out_list  = plu_mask_out_list

        return x_detector,  y_detector, xylenslight, kwargs_data_joint




    def numerics(self, supersampling_factor=4, supersampling_convolution=False):
        """
        numerical option of the image data and convolution
        :param supersampling_factor: int, factor of higher resolution sub-pixel sampling of surface brightness
        :param supersampling_convolution: bool, if True, performs (part of) the convolution on the super-sampled grid/pixels
        :return: numerical arguments
        """
        kwargs_numerics={'supersampling_factor': supersampling_factor, 'supersampling_convolution': supersampling_convolution}
        return kwargs_numerics





    def plot_segmentation(self,image_data,segments_deblend,xcenter,ycenter,c_index,font_size=20):
        """
        show segmentation map of image_data
        :param image_data:
        :param segments_deblend:
        :param nobjs:
        :param xcenter:
        :param ycenter:
        :param c_index:
        :return:
        """
        n_max=1 #np.max(image_data)*0.5
        n_min = -2#np.min(image_data)
        fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(10,6))
        ax1.imshow(np.log10(image_data), origin='lower',cmap="gist_heat")#, vmax=n_max,vmin=n_min)
        ax1.set_title('Input Image',fontsize =font_size )

        ax2.imshow(segments_deblend, origin='lower')
        for i in range(len(xcenter)):
            ax2.text(xcenter[i], ycenter[i], repr(i), color='r',size = 12)
        ax2.text(image_data.shape[0]*0.5,image_data.shape[0]*0.1,'Seg '+repr(c_index)+' '+'in center',size=12,color='white')
        ax2.set_title('Segmentations (S/N >'+repr(self.snr)+')',fontsize =font_size)


        plt.show()
        return 0



    def plot_prodata_psf(self,font_size=28,img_name='prodata_psf.pdf',img_id=0):
        """
        plot original data, processed one (data+lens light+pollution), lens light, mask, psf
        :return:
        """
        rawimage = self.raw_image
        dataimage = self.data
        len_mask = self.lens_mask
        plu_mask_out = self.plu_mask

        fig, (ax1, ax2, ax3, ax4,ax5) = plt.subplots(1, 5, figsize=(19, 10))
        ax1.imshow((rawimage), origin='lower', cmap="gist_heat")
        ax1.set_title('Original Image', fontsize=font_size)
        ax1.text(rawimage.shape[0] * 0.55, rawimage.shape[0] * 0.8, 'ID='+repr(img_id), size=12, color='white',
                 weight="bold")
        ax1.text(rawimage.shape[0] * 0.2, rawimage.shape[0] * 0.05, 'observation', size=20, color='white', weight="bold")
        ax1.axis('off')
        #
        ax2.imshow((dataimage), origin='lower', cmap="gist_heat")
        ax2.set_title('Image Data', fontsize=font_size)
        ax2.text(dataimage.shape[0] * 0.2, dataimage.shape[0] * 0.05, 'image data', size=20, color='white', weight="bold")
        ax2.axis('off')
        #
        ax3.imshow(len_mask, origin='lower')
        ax3.set_title('Lens light', fontsize=font_size)
        ax3.axis('off')
        #
        ax4.imshow(plu_mask_out, origin='lower')
        ax4.set_title('Mask', fontsize=font_size)
        ax4.axis('off')
#
        psf=self.psf
        ax5.imshow(np.log10(psf), origin='lower', cmap="gist_heat")
        ax5.set_title('lg(PSF)', fontsize=font_size)
        ax5.axis('off')

        plt.show()
        fig.savefig(img_name)
        return 0





    def plot_data_assemble(self,kwargs_seg, add_mask ,img_name='data.pdf',cutout_text='lensed image',font_size=28):
        """
        plot data assemle procession
        :param kwargs_seg: selected segmentation kwargs
        :param add_mask: int, numbers of pixels added surround the selected segmentation
        :param img_name: string, figure name
        :param cutout_text: string, figure label
        :param font_size: int, title font size
        :return:
        """
        mask = self.data_mask
        image = self.raw_image
        picked_data = self.data
        selem = np.ones((add_mask, add_mask))
        img_mask = ndimage.binary_dilation(mask.astype(np.bool), selem)
        fig, (ax1, ax2, ax3,ax4) = plt.subplots(1, 4, figsize=(19, 10))
        ax1.imshow(image, origin='lower', cmap="gist_heat")
        ax1.set_title('Cutout Image',fontsize =font_size)
        ax1.text(image.shape[0] * 0.2, image.shape[0] * 0.05, cutout_text,size=20, color='white',weight="bold")
        ax1.axis('off')
        segments_deblend_list, xcenter, ycenter, c_index=kwargs_seg
        ax2.imshow(segments_deblend_list, origin='lower')
        for i in range(len(xcenter)):
            ax2.text(xcenter[i] * 1.1, ycenter[i], 'Seg' + repr(i),  size=20,color='w',weight="bold")
        ax2.text(image.shape[0] * 0.2, image.shape[0] * 0.9, 'Seg' + repr(c_index) + ' ' + 'in center',
                 size=20, color='white',weight="bold")
        ax2.set_title('Segmentations',fontsize =font_size)
        ax2.axis('off')
        ax3.imshow(img_mask+mask, origin='lower',cmap="gist_heat")
        ax3.set_title('Selected pixels',fontsize =font_size)
        ax3.text(image.shape[0] * 0.1, image.shape[0] * 0.05, 'pixels (S/N >' + repr(self.snr) + ')',size=20, color='white',weight="bold")
        ax3.text(image.shape[0] * 0.1, image.shape[0] * 0.9, 'additional pixels', size=20, color='r',weight="bold")
        ax3.axis('off')
        ax4.imshow(picked_data, origin='lower',cmap="gist_heat")
        ax4.set_title('Processed Image',fontsize =font_size)
        ax4.axis('off')
        plt.show()
        fig.savefig(img_name)
        return 0


