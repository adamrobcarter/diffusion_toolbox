import numpy as np
import os

def edf_read_dirty(filename, headerlength=512, imsize=[2048, 2048], \
    roi=[0, 0, 2048, 2048], datatype=np.uint32, accfac=1):
    '''
    dirty and nasty edf_read for short term use!
    simple_edf_read(filename, headerlength=512, imsize=[2048,2048], roi=[0, 750, 2048, 1], datatype=np.uint32)
    imsize = xsize, ysize
    roi is x start, y start, xrange, yrange
    accfac = 1 - if this is an integer > 1, divide by this to convert to 16bit
    '''
    
    # open the file
    f = open(filename, 'rb')
    # how long is one value?
    bitlength = np.dtype(datatype).itemsize
    # skip to the start
    f.seek(headerlength+(bitlength*imsize[0]*roi[1]))
    # read the block (full lines)
    #a = f.read(bitlength*imsize[0]*roi[3])
    # sort into an image
    #b = np.fromstring(a, datatype)
    #c = np.reshape(b, (800, 2048))
    # as one line
    c = np.reshape(np.fromstring(f.read(bitlength*imsize[0]*roi[3]), datatype), (roi[3], imsize[0]))
    # close the file
    f.close()
    # crop this for the x roi
    c = c[:, roi[0]:(roi[0]+roi[2])]
    # convert to 16 bit ?
    if (datatype==np.uint32) and (accfac>1):
        c = np.uint16(c/accfac)
    
    return c
    

def vol_read_dirty_OLD(filename, volsize=None, roi=None, datatype=None):
    '''
    vol_read_dirty(filename, volsize=[2048, 2048, 1], roi=[0,0,0,2048,2048,1], datatype=np.float32)
    use the same data reordering as Henry Proudhon
    For ROIs, should probably just use a np memmap and let np deal with everything
    '''
    # open the file
    f = open(filename, 'rb')
    # if no volsize given, try to read the .info
    if volsize==None:
        volsize, datatype = read_volsize(filename)
            
    # if no roi given, try to read first slice only 
    if roi==None:
        print("reading whole volume - if you want less, specify the roi")
        roi = [0,0,0,volsize[0],volsize[1],volsize[2]]
    elif roi=='centre':
        print("read the centre slice in z")
        roi = [0,0,int(np.floor(volsize[2]/2)),volsize[0],volsize[1],1]
        
            
    # if no datatype is given, infer from the filename
    if datatype==None:
        if filename[-4::]=='.vol':
            datatype = np.float32
        elif filename[-4::]=='.raw':
            datatype = np.uint8
        else:
            datatype = np.float32
            print("guessing datatype is float")
    
    # how long is one value?
    bitlength = np.dtype(datatype).itemsize
    # skip n slices to the starting point
    f.seek(bitlength*volsize[0]*volsize[1]*roi[2])
    # read the block of m whole slices
    c = np.reshape(np.fromstring(f.read(bitlength*volsize[0]*volsize[1]*roi[5]), datatype), (roi[5], volsize[1], volsize[0]))
    # close the file
    f.close()
    # reorder the data
    c = c.transpose(2,1,0)
    # crop to the x,y roi
    c = c[roi[0]:(roi[0]+roi[3]), roi[1]:(roi[1]+roi[4]), :]
    
    return c


def vol_read_dirty(filename, volsize=None, roi=None, datatype=None):
    '''
    vol_read_dirty(filename, volsize=[2048, 2048, 1], roi=[0,0,0,2048,2048,1], datatype=np.float32)
    use the same data reordering as Henry Proudhon
    For ROIs, should probably just use a np memmap and let np deal with everything
    '''

    # if no volsize given, try to read the .info
    if volsize==None:
        volsize, datatype = read_volsize(filename)
            
    # if no roi given, try to read first slice only 
    if roi=='centre':
        print("read the centre slice in z")
        roi = [0,0,int(np.floor(volsize[2]/2)),volsize[0],volsize[1],1]
    
    # get a memmap
    mm = vol_read_virtual(filename)
    
    # read
    if roi==None:
        print("returning the whole volume")
        vol = mm[:] * 1
    else:
        vol = mm[roi[0]:(roi[0]+roi[3]), roi[1]:(roi[1]+roi[4]), roi[2]:(roi[2]+roi[5])] * 1

    return vol

def vol_read_virtual(filename, mode='r'):
    '''
    Read as a numpy memory map
    '''
    # if no volsize given, try to read the .info
    volsize, datatype = read_volsize(filename)
    mm = np.memmap(filename, datatype, mode, 0, tuple(volsize), 'F')
    return mm
            
def vol_write(data, filename):
    '''
    Write, using a numpy memory map
    '''
    if data.dtype == np.float64:
        print('64bit data... converting to float32!')
        data = data.astype(np.float32)
    if data.dtype not in [np.float32, np.uint16, np.uint8]:
        print('Strange datatype : %s - are you sure this is right?' % str(data.dtype))
    else:
        print('Datatype : %s - seems ok' % str(data.dtype))
    volsize = data.shape
    # open memory map
    mm = np.memmap(filename, data.dtype, 'write', 0, tuple(volsize), 'F')
    mm[:] = data[:]
    mm.flush()
    del mm
    print('file written')
    f = open(filename + '.info', 'wt')
    f.write('! PyHST VOLUME INFO FILE\n')
    f.write('NUM_X = %d\n' % volsize[0])
    f.write('NUM_Y = %d\n' % volsize[1])
    f.write('NUM_Z = %d\n' % volsize[2])
    f.close()
    print('.info file written')
   
def read_volsize(filename):
    if True: # try:
        pars = {}
        if filename[-3:] in ['vol', 'raw']:
            filename = filename + '.info'
        else:
            print("strange filename...  trying anyway")
            filename = filename + '.info'
        f2 = open(filename, 'rt')
        lines = f2.readlines()
        f2.close()
        for line in lines:
            parts = line.split("=")
            if len(parts)==2:
                name = parts[0].strip()
                try:
                    val = float(parts[1].strip())
                    pars.update({name:val})
                except:
                    #print("failed to read %s" % name)
                    # read as string
                    pars.update({name:parts[1]})
        if 'NUM_Z' in pars.keys():
            volsize = [pars['NUM_X'], pars['NUM_Y'], pars['NUM_Z']]
        else:
            volsize = [pars['NUM_X'], pars['NUM_Y'], 1]
        volsize = [int(x) for x in volsize]
        # calculate the data format
        fsize = os.path.getsize(filename[0:-5])
        nbytes = fsize / (volsize[0] * volsize[1] * volsize[2])
        if nbytes==1:
            datatype = np.uint8
        elif nbytes==2:
            datatype = np.uint16
        elif nbytes==4:
            datatype = np.float32
        else:
            datatype = None
            print("File size does not make sense!")
    #except:
    #    print("Failed to read .info file")
    #    volsize = None
    #    datatype = None
    return volsize, datatype

