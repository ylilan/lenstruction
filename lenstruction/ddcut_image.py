import numpy as np
import matplotlib.pyplot as plt


def sub_bkg(img, plot=False):
    from astropy.stats import SigmaClip
    from photutils import Background2D, SExtractorBackground
    sigma_clip = SigmaClip(sigma=3., iters=10)
    bkg_estimator = SExtractorBackground()
    from photutils import make_source_mask
    mask_0 = make_source_mask(img, snr=2, npixels=5, dilate_size=11)
    mask_1 = (np.isnan(img))
    mask = mask_0 + mask_1
    bkg = Background2D(img, (50, 50), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
                       mask=mask)
    from matplotlib.colors import LogNorm
    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_subplot(1, 1, 1)
    ax.imshow(img, norm=LogNorm(), origin='lower')
    # bkg.plot_meshes(outlines=True, color='#1f77b4')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    if plot:
        plt.show()
    else:
        plt.close()
    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_subplot(1, 1, 1)
    ax.imshow(mask, origin='lower')
    # bkg.plot_meshes(outlines=True, color='#1f77b4')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    if plot:
        plt.show()
    else:
        plt.close()
    back = bkg.background * ~mask_1
    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_subplot(1, 1, 1)
    ax.imshow(back, origin='lower', cmap='Greys_r')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    if plot:
        plt.show()
    else:
        plt.close()
    return img - back, back




