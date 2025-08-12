
import os
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import pandas as pd
import itk
import skimage


from sklearn.linear_model import LinearRegression
from scipy import signal
from scipy import interpolate

default_anatname='SAGT13DMPRAGE.nii'

def alignanatomical(taskname, anatname):

    rawfile, anatfile, anatalignedfile, finalfile = getimgnames(taskname, anatname)

    targetres = 200

    funcimg = nib.load(rawfile)
    functransform = funcimg.affine
    funcsize = funcimg.header.get_data_shape()
    zs = funcsize[2]
    xs = funcsize[0]
    ys = funcsize[1]

    anatimg = nib.load(anatfile)
    anattransform = anatimg.affine
    anatdata = anatimg.get_fdata()
    anatinv = np.linalg.inv(anattransform)
    anatsize = anatimg.header.get_data_shape()
    xas = anatsize[0]
    yas = anatsize[1]
    zas = anatsize[2]
    
    # Desired coordinates in native space
    print("Getting Coordinates in Native Space")
    totpts = targetres*targetres*zs
    tstarr = np.zeros((targetres,targetres,zs))
    for i in range(targetres):
        tstarr[i,:,:] = float(i)/float(targetres) * float(xs)
    xpts = np.reshape(tstarr,-1,order='F')
    tstarr.fill(0)
    for j in range(targetres):
        tstarr[:,j,:] = float(j)/float(targetres) * float(ys)
    ypts = np.reshape(tstarr,-1,order='F')
    tstarr.fill(0)
    for k in range(zs):
        tstarr[:,:,k] = k
    zpts = np.reshape(tstarr,-1,order='F')

    
    # Transform to scanner space
    print("Transforming to Scanner Space")
    
    xcoords = functransform[0,0]*xpts + functransform[0,1]*ypts + \
              functransform[0,2]*zpts + functransform[0,3]
    ycoords = functransform[1,0]*xpts + functransform[1,1]*ypts + \
              functransform[1,2]*zpts + functransform[1,3]
    zcoords = functransform[2,0]*xpts + functransform[2,1]*ypts + \
              functransform[2,2]*zpts + functransform[2,3]

    # Transform to anatomical space
    print("Transforming to Anatomical Space")
    xapts = anatinv[0,0]*xcoords + anatinv[0,1]*ycoords + \
            anatinv[0,2]*zcoords + anatinv[0,3]
    yapts = anatinv[1,0]*xcoords + anatinv[1,1]*ycoords + \
            anatinv[1,2]*zcoords + anatinv[1,3]
    zapts = anatinv[2,0]*xcoords + anatinv[2,1]*ycoords + \
            anatinv[2,2]*zcoords + anatinv[2,3]
    
    # Interpolate
    print("Interpolating")
    xgpts = np.linspace(0,xas-1,xas)
    ygpts = np.linspace(0,yas-1,yas)
    zgpts = np.linspace(0,zas-1,zas)
    gpts = (xgpts,ygpts,zgpts)
    rpts = np.zeros((totpts,3))
    rpts[:,0] = xapts
    rpts[:,1] = yapts
    rpts[:,2] = zapts
    result = interpolate.interpn(gpts,anatdata,rpts,bounds_error=False,fill_value=0)
    result = np.reshape(result,(zs,targetres,targetres))
    result = np.transpose(result,(2,1,0))

    # breakpoint()
                        
    # Save (fix resolution in affine transformation)
    print("Saving")
    affine = functransform
    affine[0,0] = affine[0,0] * float(xs)/float(targetres)
    affine[1,1] = affine[1,1] * float(ys)/float(targetres)
    tosave = nib.Nifti1Image(result,affine)
    nib.save(tosave,anatalignedfile)

    # Done
    return
   

def createimage(taskname, anatname):

    rawfile, anatfile, anatalignedfile, finalfile = getimgnames(taskname, anatname)

    # Makes final image


    threshl = 6.0
    threshh = 16.0

    funcimg = nib.load(rawfile)
    funcsize = funcimg.header.get_data_shape()
    funcdata = funcimg.get_fdata()
    funcdata = np.transpose(funcdata, (1,0,2))
    funcdata = np.flip(funcdata, (0,1))
    xs = funcsize[1]
    ys = funcsize[0]

    anatimg = nib.load(anatalignedfile)
    anatsize = anatimg.header.get_data_shape()
    anatdata = anatimg.get_fdata()
    anatdata = np.transpose(anatdata, (1,0,2))
    anatdata = np.flip(anatdata, (0,1))

 
    xas = anatsize[1]
    yas = anatsize[0]
    zs = anatsize[2]

    # Expand funcimg to size of anat
    efuncdata = np.zeros((xas,yas,zs))
    for k in range(zs):
        efuncdata[:,:,k] = skimage.transform.resize(funcdata[:,:,k],(xas,yas),order=0)

    # Find mosaic dimensions
    ximgnum = int(zs/np.sqrt(zs))
    yimgnum = int(zs/ximgnum)
    if (ximgnum*yimgnum < zs):
        yimgnum += 1

    ximgnum = int(ximgnum)
    yimgnum = int(yimgnum)

    # Scale anatomical (grayscale) and make image in mosaic
    anatdata = anatdata / np.max(anatdata) * 255
    anatimage = np.zeros((ximgnum*xas,yimgnum*yas,3))
    for k in range(zs):
        # xstart = int((k % ximgnum) * xas)
        # ystart = int(int(k / ximgnum) * yas)
        ystart = int((k % yimgnum) * yas)
        xstart = int(int(k / yimgnum) * xas)
        anatimage[xstart:xstart+xas,ystart:ystart+yas,0] = anatdata[:,:,k]
        anatimage[xstart:xstart+xas,ystart:ystart+yas,1] = anatdata[:,:,k]
        anatimage[xstart:xstart+xas,ystart:ystart+yas,2] = anatdata[:,:,k]
        
    
    # Scale functional (yellow-to-red colorscale and make image in mosaic

    funcimage = np.zeros((ximgnum*xas,yimgnum*yas,3))
    for k in range(zs):
         # xstart = int((k % ximgnum) * xas)
         # ystart = int(int(k / ximgnum) * yas)
         ystart = int((k % yimgnum) * yas)
         xstart = int(int(k / yimgnum) * xas)
         funcslice = efuncdata[:,:,k]
         if np.sum(funcslice) > 0: # Make sure there are some nonzero points in the slice
             # Red is 255 everywhere there are nonzero points
             redslice = np.zeros((xas,yas))
             redslice[funcslice != 0] = 255
             # Green is reverse-scaled from 0 to 255
             greenslice = np.zeros((xas,yas))
             funcpts = funcslice[funcslice != 0]
             funcpts[funcpts > threshh] = threshh
             funcpts = funcpts - threshl
             funcpts = funcpts / (threshh - threshl) * 255
             funcpts = 255 - funcpts
             greenslice[funcslice != 0] = funcpts
             # Blue is 0 everywhere
             funcimage[xstart:xstart+xas,ystart:ystart+yas,0] = redslice
             funcimage[xstart:xstart+xas,ystart:ystart+yas,1] = greenslice
          

    # Save anatomy and function for debugging
    # skimage.io.imsave("Functestimg.tiff",np.asarray(funcimage,'int8'))
    # skimage.io.imsave("Anattestimg.tiff",np.asarray(anatimage,'int8'))
    


    # Make overlay image
    overlayimage = anatimage
    overlayimage[:,:,0] = np.where(funcimage[:,:,0] == 0,anatimage[:,:,0],funcimage[:,:,0])
    overlayimage[:,:,1] = np.where(funcimage[:,:,0] == 0,anatimage[:,:,1],funcimage[:,:,1])
    overlayimage[:,:,2] = np.where(funcimage[:,:,0] == 0,anatimage[:,:,2],funcimage[:,:,2])

    
    # Save overlay
    skimage.io.imsave(finalfile,np.asarray(overlayimage,'int8'))

    #Done
    return
    

def thresholdcluster(taskname, anatname):

    # Makes filtered, clustered image ready to overlay on anatomical

    origfile, refnumfile, alignedfile, paramfile, costfile, strippedfile, Zfile = getfilenames(taskname)

    rawfile, anatfile, anatalignedfile, finalfile = getimgnames(taskname, anatname)

    # Thresholds
    Zthresh = 6.0
    clust = 100

    # Load in image
    img = nib.load(Zfile)
    data = img.get_fdata()
    mask = data != 0
    mask = mask.astype(float)
    header = img.header
    size = header.get_data_shape()
    xs = size[0]
    ys = size[1]
    zs = size[2]

    # Make Gaussian kernel, sigma = 1 voxel
    kern = np.zeros((5,5,5))
    for k in range(5):
        for j in range(5):
            for i in range(5):
                x = i-2
                y = j-2
                z = k-2
                kern[i,j,k] = -(x**2 + y**2 + z**2)/2
    kern = np.exp(kern)
 

    # Convolve with rms normalization, only values where originally values (no "bleeding")
    cdata = signal.convolve(data, kern, mode='same')    
    norm = np.sqrt(signal.convolve(mask, kern**2, mode='same'))
    fdata = np.zeros((xs,ys,zs))
    fdata[mask == 1] = cdata[mask==1]/norm[mask==1]


    # Threshold and cluster
    actmask = fdata >= Zthresh
    actmask = actmask.astype(float)
    stilltodo = actmask
    clusters = np.zeros((xs,ys,zs))
    while np.sum(stilltodo) > 0:
        pts = np.nonzero(stilltodo)
        xpts = pts[0]
        ypts = pts[1]
        zpts = pts[2]
        xpt = xpts[0]
        ypt = ypts[0]
        zpt = zpts[0]
        tfill = skimage.segmentation.flood_fill(actmask,(xpt,ypt,zpt),2)
        clusters[tfill == 2] = np.sum(tfill == 2)
        stilltodo[tfill == 2] = 0
    clustmask = clusters >= clust
    clustmask = clustmask.astype(float)
    result = np.multiply(fdata,clustmask)

    # Save Clusters for Debugging
    # tosave = nib.Nifti1Image(clusters,img.affine)
    # nib.save(tosave,"Clusters.nii")

    # Save Result
    tosave = nib.Nifti1Image(result,img.affine)
    nib.save(tosave,rawfile)

    # alldone
    return

    
    

def processstandard(taskname,timename,Zname=''):


    origfile, refnumfile, alignedfile, paramfile, costfile, strippedfile, Zfile = getfilenames(taskname)
    if Zname != '':
        Zfile = Zname+'Z.nii'


    # loading in data, time, and motion
    img = nib.load(strippedfile)
    data = img.get_fdata()
    affine = img.affine
    timeseries = pd.read_csv(timename,header=None)
    timec = timeseries.to_numpy()
    timec = timec.transpose()
    timec = timec[0,:]
    motparms = np.loadtxt(paramfile)

    # binarymask
    maskthresh = 0.2
    timesum = np.sum(data,3)
    maxvalue = np.max(timesum)
    mask = np.divide(timesum,maxvalue)
    
    # getting necessary information and making result array
    header = img.header
    size = header.get_data_shape()
    xs = size[0]
    ys = size[1]
    zs = size[2]
    ts = size[3]
    result = np.zeros((xs,ys,zs))

    # making regression array
    regr = np.zeros((3+12,ts))
    for i in range(ts):
        regr[0,i] = timec[i]
        regr[1,i] = i
        regr[2,i] = i**2
    regr[3:,:] = motparms
    regr = regr.transpose()
    #np.savetxt("TstRegression.txt",regr,delimiter='\t')

    # here we go
    for k in range(zs):
        print(k)
        for j in range(ys):
            for i in range(xs):
                if mask[i,j,k] > maskthresh:
                    line = data[i,j,k,:]
                    line = line.reshape(-1)
                    good = (line > 0) & (regr[:,0] >= 0)
                    regr1 = regr[good,:]
                    line1 = line[good]
                    rgrsh = regr1.shape
                    reg = LinearRegression().fit(regr1,line1)
                    resid = line1 - reg.predict(regr1)
                    residerr2 = np.sum(resid**2) / (rgrsh[0] - rgrsh[1] - 1)
                    covmat = np.linalg.inv(regr.T @ regr) * residerr2
                    stderr = np.sqrt(covmat[0,0])
                    coeff = reg.coef_
                    T = coeff[0]/stderr
                    result[i,j,k] = T
                

    # save results
    tosave = nib.Nifti1Image(result,affine)
    nib.save(tosave,Zfile)

    # alldone
    return


def processhush(taskname,timename,hushfac):

    # Still to do:
    # Sanity checking each hush factor
    origfile, refnumfile, alignedfile, paramfile, costfile, strippedfile, Zfile = getfilenames(taskname)

    # loading in data and text file
    img = nib.load(strippedfile)
    data = img.get_fdata()
    affine = img.affine
    timeseries = pd.read_csv(timename,header=None)
    timec = timeseries.to_numpy()
    timec = timec.transpose()
    timec = timec[0,:]
    motparms = np.loadtxt(paramfile)


    # binarymask
    maskthresh = 0.2
    timesum = np.sum(data,3)
    maxvalue = np.max(timesum)
    mask = np.divide(timesum,maxvalue)
    
    # getting necessary information and making result array
    header = img.header
    size = header.get_data_shape()
    xs = size[0]
    ys = size[1]
    zs = size[2]
    ts = size[3]
    result = np.zeros((xs,ys,zs))

    # making regression array
    repnum = int(ts/hushfac)
    regr = np.zeros((hushfac,repnum,int(1+3*hushfac)))
    for i in range(repnum):
        for j in range(hushfac):
            regr[:,i,0] = timec[i] # Contrast
            regr[j,i,j+1] = 1 # Constant
            regr[j,i,j+1+hushfac] = i # Linear
            regr[j,i,j+1+2*hushfac] = i**2 # Quadratic
    regr = np.transpose(regr,(1,0,2))
    regr = np.reshape(regr,(int(hushfac*repnum),int(1+3*hushfac)))
    # insert motion
    regr1 = np.zeros((ts,int(1+3*hushfac+12)))
    regr1[:,:1+3*hushfac] = regr
    regr1[:,1+3*hushfac:] = np.transpose(motparms)
    regr = regr1

    # Print out regression array to debug
    # np.savetxt('HushRegression.txt',regr,delimiter='\t')
    # return
 
    # here we go
    for k in range(zs):
        print(k)
        for j in range(ys):
            for i in range(xs):
                if mask[i,j,k] > maskthresh:
                    line = data[i,j,k,:]
                    line = line.reshape(-1)
                    good = (line > 0) & (regr[:,0] >= 0)
                    regr1 = regr[good,:]
                    line1 = np.log(line[good])
                    rgrsh = regr1.shape
                    reg = LinearRegression(fit_intercept=False).fit(regr1,line1)
                    resid = line1 - reg.predict(regr1)
                    residerr2 = np.sum(resid**2) / (rgrsh[0] - rgrsh[1])
                    covmat = np.linalg.inv(regr.T @ regr) * residerr2
                    stderr = np.sqrt(covmat[0,0])
                    coeff = reg.coef_
                    T = coeff[0]/stderr
                    result[i,j,k] = T
                

    # save results
    tosave = nib.Nifti1Image(result,affine)
    nib.save(tosave,Zfile)

    # alldone
    return

def alignimage(fixedimg,movingimg):

    # Align movingimg to fixedimg, this is just for testing

    # Get in itk format
    print("In routine, putting into itk")
    
    fixed_itk = itk.image_from_array(fixedimg)
    moving_itk = itk.image_from_array(movingimg)

    # Set stuff up
    print("Setting stuff up")
    par_obj = itk.ParameterObject.New()
    par_map = par_obj.GetDefaultParameterMap('affine')
    par_obj.AddParameterMap(par_map)

    # Call the routine
    print("Calling routine")
    result_itk, rtp = itk.elastix_registration_method(fixed_itk, moving_itk,
                                    parameter_object=par_obj)


    # Return the result
    print("Done, back into numpy")
    resultimg = itk.array_from_image(result_itk)
    return resultimg, rtp


def alignfunc(taskname):

    origfile, refnumfile, alignedfile, paramfile, costfile, strippedfile, Zfile = getfilenames(taskname)

    # Loading data
    print("Loading data")
    img = nib.load(origfile)
    data = img.get_fdata()
    header = img.header
    size = header.get_data_shape()
    ts = size[3]

    # Loading refnum
    with open(refnumfile, "r") as file:
        refnum = file.readline()
    refnum = int(refnum)

    # Setting up
    print("Setting up")
    par_obj = itk.ParameterObject.New()
    par_map = par_obj.GetDefaultParameterMap('affine')
    par_map["ResampleInterpolator"] = ('FinalLinearInterpolator',)
    par_obj.AddParameterMap(par_map)
    fitparams = np.zeros((12,ts))
    fixed_itk = itk.image_from_array(data[:,:,:,refnum])

    # Running the routine
    for i in range(ts):
        if i != refnum:
            print('Aligning %d to %d' % (i,refnum))
            moving_itk = itk.image_from_array(data[:,:,:,i])
            result_itk, rtp = itk.elastix_registration_method(fixed_itk, moving_itk,
                                    parameter_object=par_obj)
            resultimg = itk.array_from_image(result_itk)
            data[:,:,:,i] = np.transpose(resultimg,(2,1,0))
            params = rtp.GetParameter(0,"TransformParameters")
            for j in range(12):
                fitparams[j,i] = float(params[j])

    # Setting identity transform for refnum
    fitparams[0,refnum] = 1.0
    fitparams[4,refnum] = 1.0
    fitparams[8,refnum] = 1.0

    print("Saving")
    tosave = nib.Nifti1Image(data,img.affine)
    nib.save(tosave,alignedfile)
    np.savetxt(paramfile,fitparams,delimiter='\t')

    return


def alignhush(taskname,hushfac):

    origfile, refnumfile, alignedfile, paramfile, costfile, strippedfile, Zfile = getfilenames(taskname)

    # Loading data
    print("Loading data")
    img = nib.load(origfile)
    data = img.get_fdata()
    header = img.header
    size = header.get_data_shape()
    ts = size[3]
    xs = size[0]
    ys = size[1]
    zs = size[2]
    ntrials = int(ts/hushfac)
    data = np.reshape(data,(xs,ys,zs,ntrials,hushfac))

    # Loading refnum
    with open(refnumfile, "r") as file:
        refnum = file.readline()
    refnum = int(refnum)

    # Setting up
    print("Setting up")
    par_obj = itk.ParameterObject.New()
    par_map = par_obj.GetDefaultParameterMap('affine')
    par_map["ResampleInterpolator"] = ('FinalLinearInterpolator',)
    par_obj.AddParameterMap(par_map)
    fitparams = np.zeros((12,ntrials,hushfac))

    # Running the routine
    for i in range(ntrials):
        if i != refnum:
            for k in range(hushfac):
                print('Aligning %d to %d' % (i*hushfac+k,refnum*hushfac+k))
                fixed_itk = itk.image_from_array(data[:,:,:,refnum,k])
                moving_itk = itk.image_from_array(data[:,:,:,i,k])
                result_itk, rtp = itk.elastix_registration_method(fixed_itk, moving_itk,
                                        parameter_object=par_obj)
                resultimg = itk.array_from_image(result_itk)
                data[:,:,:,i,k] = np.transpose(resultimg,(2,1,0))
                params = rtp.GetParameter(0,"TransformParameters")
                for j in range(12):
                    fitparams[j,i,k] = float(params[j])

    # Setting identity transform for refnum
    fitparams[0,refnum,:] = 1.0
    fitparams[4,refnum,:] = 1.0
    fitparams[8,refnum,:] = 1.0

    print("Saving")
    data = np.reshape(data,(xs,ys,zs,ts))
    fitparams = np.reshape(fitparams,(12,ts))
    tosave = nib.Nifti1Image(data,img.affine)
    nib.save(tosave,alignedfile)
    np.savetxt(paramfile,fitparams,delimiter='\t')

    return



def stripbadframes(taskname,thresh):

    origfile, refnumfile, alignedfile, paramfile, costfile, strippedfile, Zfile = getfilenames(taskname)

    # Strips motion-excessive frames, writing "0"

    # Load data and get info
    img = nib.load(alignedfile)
    data = img.get_fdata()
    header = img.header
    size = header.get_data_shape()
    costs = np.loadtxt(costfile)
    size = header.get_data_shape()

    for i in range(size[3]):
        if costs[i] > thresh:
            data[:,:,:,i] = 0

    tosave = nib.Nifti1Image(data,img.affine)
    nib.save(tosave,strippedfile)
    
    
   
def getfinalcost(taskname):

    origfile, refnumfile, alignedfile, paramfile, costfile, strippedfile, Zfile = getfilenames(taskname)
    
    # Load data and get info
    img = nib.load(alignedfile)
    data = img.get_fdata()
    header = img.header
    size = header.get_data_shape()

    # Loading refnum
    with open(refnumfile, "r") as file:
        refnum = file.readline()
    refnum = int(refnum)


    ts = size[3]
    cost = np.zeros(ts)

    # Compute costs based on reference
    divis = np.sum(data[:,:,:,refnum]**2)
    for i in range(ts):
        cost[i] = np.sum(np.subtract(data[:,:,:,i],data[:,:,:,refnum])**2) / divis

    # Write out results
    np.savetxt(costfile,cost)
    
    return


def getfinalcosthush(taskname,hushfac):

    origfile, refnumfile, alignedfile, paramfile, costfile, strippedfile, Zfile = getfilenames(taskname)
    
    # Load data and get info
    img = nib.load(alignedfile)
    data = img.get_fdata()
    header = img.header
    size = header.get_data_shape()
    xs = size[0]
    ys = size[1]
    zs = size[2]

    # Loading refnum
    with open(refnumfile, "r") as file:
        refnum = file.readline()
    refnum = int(refnum)


    ts = size[3]
    ntrials = int(ts/hushfac)
    
    cost = np.zeros((ntrials,hushfac))
    data = np.reshape(data,(xs,ys,zs,ntrials,hushfac))
    

    # Compute costs based on reference
    divis = np.zeros(hushfac)
    for i in range(hushfac):
        divis[i] = np.sum(data[:,:,:,refnum,i]**2)
        
    for i in range(ntrials):
        for j in range(hushfac):
            cost[i,j] = np.sum(np.subtract(data[:,:,:,i,j],data[:,:,:,refnum,j])**2) / divis[j]

    # Write out results
    cost = np.reshape(cost,ts)
    np.savetxt(costfile,cost)
    
    return


def getfilenames(taskname):

    # Getting all the filenames for a given task
    origfile = taskname + ".nii"
    refnumfile = taskname + "Refnum.txt"
    alignedfile = taskname + "Aligned.nii"
    paramfile = taskname + "fitparams.txt"
    costfile = taskname + "FinalCost.txt"
    strippedfile = taskname + "Stripped.nii"
    Zfile = taskname + "Z.nii"

    return origfile, refnumfile, alignedfile, paramfile, costfile, strippedfile, Zfile
      

def getimgnames(taskname, anatname):

    # Getting all the imagenames

    rawfile = taskname + "FuncResults.nii"
    anatfile = anatname
    finalfile = taskname + "OverlayResults.tiff"
    anatalignedfile = taskname + "AlignedAnatomy.nii"

    return rawfile, anatfile, anatalignedfile, finalfile

def findreference(taskname):

    # Finds best reference frame

    origfile, refnumfile, alignedfile, paramfile, costfile, strippedfile, Zfile = getfilenames(taskname)

    # Load data and get info
    img = nib.load(origfile)
    data = img.get_fdata()
    header = img.header
    size = header.get_data_shape()
    ts = size[3]
    xs = size[0]
    ys = size[1]
    zs = size[2]

    # Compute correlations
    data = np.reshape(data,(xs*ys*zs,ts))
    data = np.transpose(data)
    diffs = np.corrcoef(data)
    diffs = 1 - diffs**2

 
    # Find optimal frame
    diff1 = np.sum(diffs,1)
    minvalue = np.argmin(diff1)
    

    # Print out
    with open(refnumfile, 'w') as file:
        file.write(str(minvalue))

    # Done
    return

def findreferencehush(taskname,hushfac):

    # Finds best reference frame

    origfile, refnumfile, alignedfile, paramfile, costfile, strippedfile, Zfile = getfilenames(taskname)

    # Load data and get info
    img = nib.load(origfile)
    data = img.get_fdata()
    header = img.header
    size = header.get_data_shape()
    ts = size[3]
    xs = size[0]
    ys = size[1]
    zs = size[2]
    ntrials = int(ts/hushfac)
    nvol = int(xs*ys*zs)

    # Compute correlations
    data = np.reshape(data,(nvol,ntrials,hushfac))
    data = np.transpose(data,(0,2,1))
    data = np.reshape(data,((nvol*hushfac),ntrials))
    data = np.transpose(data)
    diffs = np.corrcoef(data)
    diffs = 1 - diffs**2

 
    # Find optimal frame
    diff1 = np.sum(diffs,1)
    minvalue = np.argmin(diff1)
    

    # Print out
    with open(refnumfile, 'w') as file:
        file.write(str(minvalue))

    # Done
    return


def fixhush(taskname,hushfac):

    origfile = taskname + ".nii"
    hushfile = taskname + "HUSH.nii"

    img = nib.load(origfile)
    data = img.get_fdata()
    header = img.header
    size = header.get_data_shape()
    ts = size[3]
    stnum = ts % hushfac

    data = data[:,:,:,stnum:]
    
    tosave = nib.Nifti1Image(data,img.affine)
    nib.save(tosave,hushfile)

    return
    
    
def verbslogtotime(logfile,timefile):

    with open(logfile, 'r') as file:
        txt = file.readlines()

    txt = txt[6:46]

    time = np.zeros(40, dtype=int)
    for i in range(40):
        parts = txt[i].split()
        tststr = parts[3]
        if tststr[0:5] == 'Verbs':
            time[i] = 1
    
    np.savetxt(timefile,time)
    return


def processall(taskname, timename, anatname=default_anatname):

    findreference(taskname)
    alignfunc(taskname)
    getfinalcost(taskname)
    stripbadframes(taskname,0.008)
    processstandard(taskname,timename)
    thresholdcluster(taskname)
    alignanatomical(taskname, anatname)
    createimage(taskname)

    return

def processmotor(taskname, lhtimename, rhtimename, anatname=default_anatname):

    findreference(taskname)
    alignfunc(taskname)
    getfinalcost(taskname)
    stripbadframes(taskname,0.008)
    processstandard(taskname,lhtimename,'LeftHand')
    processstandard(taskname,rhtimename,'RightHand')
    thresholdcluster('LeftHand', anatname)
    thresholdcluster('RightHand', anatname)
    alignanatomical('LeftHand', anatname)
    alignanatomical('RightHand', anatname)
    createimage('LeftHand')
    createimage('RightHand')
    
    



def hushprocessall(taskname, timename, hushfac, anatname=default_anatname):

    
    fixhush(taskname,hushfac)
    hushname = taskname + "HUSH"
    findreferencehush(hushname,hushfac)
    alignhush(hushname,hushfac)
    getfinalcosthush(hushname,hushfac)
    stripbadframes(hushname,0.008)
    processhush(hushname,timename,hushfac)
    thresholdcluster(hushname, anatname)
    alignanatomical(hushname, anatname)
    createimage(hushname, anatname)

    return



