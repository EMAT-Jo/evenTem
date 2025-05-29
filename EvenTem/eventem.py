"""
copyright University of Antwerp 2025

author: Arno Annys
"""
import os
base_dir = os.path.dirname(__file__)
import sys
sys.path.append(os.path.join(base_dir,'binaries'))
import eventem
import numpy as np
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore", message="The value of the smallest subnormal for <class 'numpy.float64'> type is zero.")

def safelog(x):
    return np.log(1+1000*x/np.max(x))

class Pacbed(eventem.Pacbed):
    """
    Position avg CBED processor

    """
    def __init__(self, nx,ny,repetitions,filename):
        """
        Instantiate a Position averaged convergend beam electron diffraction (PACBED) processor

        Parameters
        ----------
        nx : int
            number of pixels in x direction
        ny : int
            number of pixels in y direction
        repetitions : int
            number of complete scans in the dataset
        filename : str
            path to the dataset file

        Returns
        -------
        None.

        """
        super().__init__(repetitions)
        self.b_cumulative = True
        self.nx = nx
        self.ny = ny    
        self.set_file(filename)

    # @property
    # def nx(self):
    #     """
    #     int : number of pixels in x direction
    #     """
    #     return super().nx
    
    # @nx.setter
    # def nx(self, value):
    #     self.nx = value

    # @property
    # def ny(self):
    #     """
    #     int : number of pixels in y direction
    #     """
    #     return super().ny
    
    # @ny.setter
    # def ny(self, value):
    #     self.ny = value
        
    
    @property
    def DwellTime(self):
        """
        int : dwell time
        """
        return super().dt
    
    @DwellTime.setter
    def DwellTime(self, value):
        self.dt = value

    @property
    def DetectorSize(self):
        """
        int : size of the detector
        """
        return super().detector_size
    
    @DetectorSize.setter
    def DetectorSize(self, value):
        self.detector_size = value


    def Run(self):
        """
        Run the virtual STEM reconstruction
        """
        super().run()

    @property
    def image(self):
        """
        2D numpy array : reconstructed image [nx,ny]
        """
        return np.array(super().Pacbed_image).reshape(self.DetectorSize,self.DetectorSize).T # transpose due to difference with C++ convention
    

    def PlotImage(self,log=False):
        """
        Plot the reconstructed image

        Parameters
        ----------
        log : bool
            plot the image in log scale
        """
        if log:
            fig, ax = plt.subplots(1,1,figsize=(10,10))
            ax.imshow(safelog(self.image),cmap='gnuplot')
        else:
            fig, ax = plt.subplots(1,1,figsize=(10,10))
            ax.imshow(self.image,cmap='gnuplot')
        ax.axis('off')

    def GetCOM(self):
        """
        Get the center of mass of the PACBED image

        Returns
        -------
        tuple : x,y coordinates of the center of mass

        """

        X,Y = np.meshgrid(np.arange(self.DetectorSize), np.arange(self.DetectorSize))
        X = X.flatten()
        Y = Y.flatten()
        I = self.image.astype(float).flatten()
        COM = np.array([np.sum(X*I)/np.sum(I), np.sum(Y*I)/np.sum(I)])
        return COM
        


class vSTEM(eventem.vSTEM):
    """
    virtual STEM processor

    """
    def __init__(self, nx,ny,repetitions,filename):
        """
        Instantiate a virtual STEM processor

        Parameters
        ----------
        nx : int
            number of pixels in x direction
        ny : int
            number of pixels in y direction
        repetitions : int
            number of complete scans in the dataset
        filename : str
            path to the dataset file

        Returns
        -------
        None.

        """
        super().__init__(repetitions)
        self.b_cumulative = True
        self.nx = nx
        self.ny = ny    
        self.set_file(filename)

    # @property
    # def nx(self):
    #     """
    #     int : number of pixels in x direction
    #     """
    #     return super().nx
    
    # @nx.setter
    # def nx(self, value):
    #     self.nx = value

    # @property
    # def ny(self):
    #     """
    #     int : number of pixels in y direction
    #     """
    #     return super().ny
    
    # @ny.setter
    # def ny(self, value):
    #     self.ny = value
        
    @property
    def DwellTime(self):
        """
        int : dwell time
        """
        return super().dt
    
    @DwellTime.setter
    def DwellTime(self, value):
        self.dt = value

    @property
    def DetectorSize(self):
        """
        int : size of the detector
        """
        return super().detector_size
    
    @DetectorSize.setter
    def DetectorSize(self, value):
        self.detector_size = value

    @property
    def InnerRadia(self):
        """
        list of int : inner radii of the annular detectors
        """
        return super().inner_radia
    
    @InnerRadia.setter
    def InnerRadia(self, value):
        self.inner_radia = value

    @property
    def OuterRadia(self):
        """
        list of int : outer radii of the annular detectors
        """
        return self.outer_radia
    
    @OuterRadia.setter
    def OuterRadia(self, value):
        self.outer_radia = value

    @property
    def Centers(self):
        """
        list of tuples : centers of the annular detectors
        """
        return self.offsets
    
    @Centers.setter
    def Centers(self, value):
        """
        list of tuples : centers of the annular detectors
        """
        self.set_offsets(value)

    def GetDetector(self):
        """
        Get the detector image

        Returns
        -------
        2D numpy array
           
        """
        return np.array(self.get_detector()).reshape(self.DetectorSize,self.DetectorSize)
    
    def PlotDetector(self,pacbed=None,log=False):
        """
        Plot the detector image
        """
        detector = self.GetDetector()
        if pacbed is None:
            fig, ax = plt.subplots(1,1,figsize=(10,10))
            ax.imshow(detector)
        elif isinstance(pacbed, Pacbed):
            if pacbed.image.shape == (self.DetectorSize,self.DetectorSize):
                self.incircle=plt.Circle(self.offsets[0],self.inner_radia[0], fill=False, color='red', linewidth=2)
                self.outcircle=plt.Circle(self.offsets[0],self.outer_radia[0], fill=False, color='red', linewidth=2)

                fig, ax = plt.subplots(1,1,figsize=(10,10))
                if log:
                    ax.imshow(safelog(pacbed.image),cmap='gray')
                else:
                    ax.imshow(pacbed.image,cmap='gray')
                ax.imshow(detector,alpha=0.2,cmap='Reds')
                ax.add_patch(self.incircle)
                ax.add_patch(self.outcircle)
                ax.axis('off')
        else:
            raise ValueError("pacbed should be a Pacbed object with the same detector size as the vSTEM object, or None")
        
    def MakeDPCDetector(self,inner_radius,outer_radius,N_segments,segment_weights,center,rotation):
        rotation_rad = rotation*np.pi/180
        _mask = np.zeros((self.DetectorSize, self.DetectorSize))
        for i in range(self.DetectorSize):
            for j in range(self.DetectorSize):
                if (i-center[0])**2 + (j-center[1])**2 < outer_radius**2:
                    if (i-center[0])**2 + (j-center[1])**2 > inner_radius**2:
                        angle = np.arctan2(i-center[0], j-center[1])
                        angle += rotation_rad
                        if angle < 0:
                            angle += 2*np.pi
                        rotation_rad = rotation*np.pi/180
                        segment = int(N_segments*angle/(2*np.pi))
                        _mask[i, j] = segment_weights[segment]
        
        _mask = _mask.flatten().astype(int)
        self.set_detector_mask(_mask)

    def Run(self):
        """
        Run the virtual STEM reconstruction
        """
        super().run()

    @property
    def image(self):
        """
        2D numpy array : reconstructed image [nx,ny]
        """
        return np.array(super().vSTEM_image).reshape(super().ny,super().nx)
    
    @property
    def image_stack(self):
        """
        3D numpy array : reconstructed image stack [nx,ny,repetitions]
        """
        return np.array(super().vSTEM_stack).reshape(super().repetitions+1,super().ny,super().nx)[:-1,:,:]
    
    def PlotImage(self):
        """
        Plot the reconstructed image
        """
        fig, ax = plt.subplots(1,1,figsize=(10,10))
        ax.imshow(self.image,cmap='gray')
        ax.axis('off')


class Var(eventem.Var):
    """
    Variance processor

    """
    def __init__(self, nx,ny,repetitions,filename):
        """
        Instantiate a variance processor

        Parameters
        ----------
        nx : int
            number of pixels in x direction
        ny : int
            number of pixels in y direction
        repetitions : int
            number of complete scans in the dataset
        filename : str
            path to the dataset file

        Returns
        -------
        None.

        """
        super().__init__(repetitions)
        self.b_cumulative = True
        self.nx = nx
        self.ny = ny    
        self.set_file(filename)

    # @property
    # def nx(self):
    #     """
    #     int : number of pixels in x direction
    #     """
    #     return super().nx
    
    # @nx.setter
    # def nx(self, value):
    #     self.nx = value

    # @property
    # def ny(self):
    #     """
    #     int : number of pixels in y direction
    #     """
    #     return super().ny
    
    # @ny.setter
    # def ny(self, value):
    #     self.ny = value
        
    @property
    def DwellTime(self):
        """
        int : dwell time
        """
        return super().dt
    
    @DwellTime.setter
    def DwellTime(self, value):
        self.dt = value

    @property
    def DetectorSize(self):
        """
        int : size of the detector
        """
        return super().detector_size
    
    @DetectorSize.setter
    def DetectorSize(self, value):
        self.detector_size = value


    @property
    def Center(self):
        """
        tuple : 
        """
        return self.offset
    
    @Center.setter
    def Center(self, value):
        """
        tuple: centers 
        """
        super().set_offset(value)

    def Run(self):
        """
        Run the virtual STEM reconstruction
        """
        super().run()

    @property
    def image(self):
        """
        2D numpy array : reconstructed image [nx,ny]
        """
        return np.array(super().Var_image).reshape(super().ny,super().nx)
    
    def PlotImage(self):
        """
        Plot the reconstructed image
        """
        fig, ax = plt.subplots(1,1,figsize=(10,10))
        ax.imshow(self.image,cmap='gray')
        ax.axis('off')
        

class Ricom(eventem.Ricom):
    """
    Ricom processor
    """
    def __init__(self, nx,ny,repetitions,filename):
        """
        Instantiate a real-time integrated center-of-mass processor

        Parameters
        ----------
        nx : int
            number of pixels in x direction
        ny : int
            number of pixels in y direction
        repetitions : int
            number of complete scans in the dataset
        filename : str
            path to the dataset file

        Returns
        -------
        None.

        """
        super().__init__(repetitions)
        self.b_cumulative = True
        self.nx = nx
        self.ny = ny    
        self.set_file(filename)
        self.n_threads = 8

    # @property
    # def nx(self):
    #     """
    #     int : number of pixels in x direction
    #     """
    #     return super().nx

    # @nx.setter
    # def nx(self, value):
    #     self.nx = value

    # @property
    # def ny(self):
    #     """
    #     int : number of pixels in y direction
    #     """
    #     return super().ny

    # @ny.setter
    # def ny(self, value):
    #     self.ny = value
        

    @property
    def DwellTime(self):
        """
        int : dwell time
        """
        return super().dt
    
    @DwellTime.setter
    def DwellTime(self, value):
        self.dt = value

    @property
    def DetectorSize(self):
        """
        int : size of the detector
        """
        return super().detector_size

    @DetectorSize.setter
    def DetectorSize(self, value):
        self.detector_size = value


    def Run(self):
        """
        Run the Ricom reconstruction
        """
        self.run()

    def SetKernel(self, kernelsize,rotation=0):
        """
        Set the kernel for the Ricom reconstruction

        Parameters
        ----------
        kernelsize : int 
            size of the kernel
        rotation : int
            rotation of the real and momentum space

        Returns
        -------
        None.

        """
        self.KS = kernelsize
        super().set_kernel(kernelsize,rotation)

    @property
    def CoMx_image(self):
        """
        2D numpy arrays : center of mass image X
        """
        return np.array(self.comx_image).reshape(self.ny,self.nx)
    
    @property
    def CoMy_image(self):
        """
        2D numpy arrays : center of mass image Y
        """
        return np.array(self.comy_image).reshape(self.ny,self.nx)
    
    @property
    def CoM(self):
        """
        tuple : center of mass
        """
        return self.offset
    
    @CoM.setter
    def CoM(self,value):
        self.set_offset(value)

    @property
    def image(self):
        """
        2D numpy array : reconstructed image [nx,ny]
        """
        return np.array(self.ricom_image).reshape(self.ny,self.nx)

    @property
    def image_stack(self):
        """
        3D numpy array : reconstructed image stack [nx,ny,repetitions]
        """
        return np.array(self.ricom_stack).reshape(self.repetitions+1,self.ny,self.nx)[:-1,:,:]

    def PlotImage(self,crop_kernel=True):
        """
        Plot the reconstructed image

        Parameters
        ----------
        crop_kernel : bool
            crop a border from the image that is the size of the kernel
        """
        if crop_kernel:
            border = self.KS
        else:
            border = 0
        fig, ax = plt.subplots(1,1,figsize=(10,10))
        ax.imshow(self.image[border:-border,border:-border])
        ax.axis('off')


class Electron(eventem.Electron):
    """
    Convert to electron file format

    """
    def __init__(self,nx,ny,repetitions,filename):
        """
        Instantiate 

        Parameters
        ----------
        nx : int
            number of pixels in x direction
        ny : int
            number of pixels in y direction
        repetitions : int
            number of complete scans in the dataset
        filename : str
            path to the dataset file

        Returns
        -------
        None.

        """
        super().__init__(repetitions)
        self.b_cumulative = True
        self.nx = nx
        self.ny = ny    
        self.set_file(filename)

    # @property
    # def nx(self):
    #     """
    #     int : number of pixels in x direction
    #     """
    #     return super().nx
    
    # @nx.setter
    # def nx(self, value):
    #     self.nx = value

    # @property
    # def ny(self):
    #     """
    #     int : number of pixels in y direction
    #     """
    #     return super().ny
    
    # @ny.setter
    # def ny(self, value):
    #     self.ny = value

    @property
    def DetectorSize(self):
        """
        int : size of the detector
        """
        return super().detector_size
    
    @DetectorSize.setter
    def DetectorSize(self, value):
        self.detector_size = value
        
    @property
    def DwellTime(self):
        """
        int : dwell time
        """
        return super().dt
    
    @DwellTime.setter
    def DwellTime(self, value):
        self.dt = value


    @property
    def xCrop(self):
        """
        int : x crop
        """
        return super().x_crop
    
    @xCrop.setter
    def xCrop(self, value):
        self.x_crop = value

    @property
    def yCrop(self):
        """
        int : y crop
        """
        return super().y_crop
    
    @yCrop.setter
    def yCrop(self, value):
        self.y_crop = value

    @property
    def ScanBin(self):
        """
        int : scan bin
        """
        return super().scan_bin
    
    @ScanBin.setter
    def ScanBin(self, value):
        self.scan_bin = value

    @property
    def DetectorBin(self):
        """
        int : detector bin
        """
        return super().detector_bin
    
    @DetectorBin.setter
    def DetectorBin(self, value):
        self.detector_bin = value

    @property
    def DropRate(self):
        """
        int : drop rate
        """
        return super().drop_rate
    
    @DropRate.setter
    def DropRate(self, value):
        self.drop_rate = value

    @property
    def ClusterRange(self):
        """
        int : cluster range
        """
        return super().cluster_range
    
    @ClusterRange.setter
    def ClusterRange(self, value):
        self.cluster_range = value

    @property
    def Dspace(self):
        """
        int : dspace
        """
        return super().dspace
    
    @Dspace.setter
    def Dspace(self, value):
        self.dspace = value

    @property
    def Dtime(self):
        """
        int : dtime
        """
        return super().dtime
    
    @Dtime.setter
    def Dtime(self, value):
        self.dtime = value

    @property
    def Decluster(self):
        """
        bool: decluster
        """
        return super().decluster
    
    @Decluster.setter
    def Decluster(self, value):
        self.decluster = value

    @property
    def Nthreads(self):
        """
        int : number of threads
        """
        return super().n_threads
    
    @Nthreads.setter
    def Nthreads(self, value):
        self.n_threads = value


    def Run(self):
        """
        Run
        """
        super().run()


class Roi(eventem.Roi):
    """
    Roi processor
    """
    def __init__(self, nx,ny,repetitions,filename,extract_4D=False):
        """
        Instantiate a ROI processor

        Parameters
        ----------
        nx : int
            number of pixels in x direction
        ny : int
            number of pixels in y direction
        repetitions : int
            number of complete scans in the dataset
        filename : str
            path to the dataset file

        Returns
        -------
        None.

        """
        super().__init__(repetitions,extract_4D)
        self.extract_4D = extract_4D
        self.b_cumulative = True
        self.nx = nx
        self.ny = ny    
        self.width = nx
        self.height = ny
        self.set_file(filename)
        self.set_roi(x=0,y=0,width=nx,height=ny) #default ROI is full image

    # @property
    # def nx(self):
    #     """
    #     int : number of pixels in x direction
    #     """
    #     return super().nx
    
    # @nx.setter
    # def nx(self, value):
    #     self.nx = value

    # @property
    # def ny(self):
    #     """
    #     int : number of pixels in y direction
    #     """
    #     return super().ny
    
    # @ny.setter
    # def ny(self, value):
    #     self.ny = value
        
    @property
    def DwellTime(self):
        """
        int : dwell time
        """
        return super().dt
    
    @DwellTime.setter
    def DwellTime(self, value):
        self.dt = value

    @property
    def ROI(self):
        """
        tuple : ROI coordinates as [x_origin,y_origin,width,height]
        """
        return np.array(super().get_roi())
    
    @ROI.setter
    def ROI(self, value):
        self.width = value[2]
        self.height = value[3]
        self.set_roi(x=value[0],y=value[1],width=value[2],height=value[3])

    def SetROI(self,x_origin,y_origin,width,height):
        """
        Set the ROI

        Parameters
        ----------
        x_origin : int
            x origin
        y_origin : int
            y origin
        width : int
            width of the ROI
        height : int
            height of the ROI

        Returns
        -------
        None.

        """
        self.width = width
        self.height = height
        self.set_roi(x=x_origin,y=y_origin,width=width,height=height)

    def SetROIMask(self,masks):
        """
        Set the ROI

        Parameters
        ----------
        x_origin : int
            x origin
        y_origin : int
            y origin
        width : int
            width of the ROI
        height : int
            height of the ROI

        Returns
        -------
        None.

        """
        self.width = self.ny #carefull with the convention here
        self.height = self.nx  #carefull with up the convention here
        mask_list = [masks[i].flatten().astype(int) for i in range(len(masks))]
        self.set_roi_mask(mask_list)

    @property
    def ScanImage(self):
        """
        2D numpy array : scan image
        """
        return np.array(super().Roi_scan_image).reshape(self.width,self.height)
    
    @property
    def DiffractionPattern(self):
        """
        2D numpy array : diffraction pattern
        """
        return np.array(super().Roi_diffraction_pattern).reshape(self.DetectorSize,self.DetectorSize)
    
    @property
    def ScanImageStack(self):
        """
        3D numpy array : scan image stack
        """
        return np.array(super().Roi_scan_image_stack).reshape(self.repetitions+1,self.width,self.height)[:-1,:,:]
    
    @property
    def DiffractionPatternStack(self):
        """
        3D numpy array : diffraction pattern stack
        """
        return np.array(super().Roi_diffraction_pattern_stack).reshape(self.repetitions+1,self.DetectorSize,self.DetectorSize)[:-1,:,:]
    
    @property
    def Roi4D(self):
        """
        4D numpy array : 4D dataset
        """
        if self.extract_4D:
            return np.array(super().get_4D(),dtype=np.uint8)
        else:
            print("4D ROI is only extracted when extract_4D is set to True")
            return None
    
    @property
    def DetectorSize(self):
        """
        int : size of the detector
        """
        return super().detector_size
    
    @DetectorSize.setter
    def DetectorSize(self, value):
        self.detector_size = value

    @property
    def DetectorBin(self):
        """
        int : detector bin
        """
        return super().det_bin
    
    @DetectorBin.setter
    def DetectorBin(self, value):
        self.det_bin = value

    def Run(self):
        """
        Run the ROI reconstruction
        """
        if self.extract_4D:
            print(f"extracting 4D sub-dataset that requires {self.width*self.height*(self.DetectorSize/self.DetectorBin)**2/1e9:.2f} GB of RAM")
        super().run()
    


class FourD():
    """
    makes a 4D dataset from events
    """

    def __init__(self,nx,ny,repetitions,filename,output_filename,bitdepth,compression_factor=7):

        if repetitions != 1:
            raise ValueError("repetitions should be 1 for FourD processor")
        
        if bitdepth not in [8,16,32]:
            raise ValueError("bitdepth should be 8, 16 or 32")
        if bitdepth == 8:
            self.super = eventem.FourD8(output_filename,repetitions,bitdepth=bitdepth,compression_factor=compression_factor)
        elif bitdepth == 16:
            self.super = eventem.FourD16(output_filename,repetitions,bitdepth=bitdepth,compression_factor=compression_factor)
        elif bitdepth == 32:
            self.super = eventem.FourD32(output_filename,repetitions,bitdepth=bitdepth,compression_factor=compression_factor)

        self.super.b_cumulative = True
        self.super.nx = nx
        self.super.ny = ny    
        self.super.set_file(filename)
        self.super.det_bin = 1 
        self.super.scan_bin = 1
        self.super.chunksize = 2


    @property
    def DetectorSize(self):
        """
        int : size of the detector
        """
        return self.super.detector_size
    
    @DetectorSize.setter
    def DetectorSize(self, value):
        self.super.detector_size = value

    @property
    def DetectorBin(self):
        """
        int : detector bin
        """
        return self.super.det_bin
    
    @DetectorBin.setter
    def DetectorBin(self, value):
        self.super.det_bin = value

    @property
    def ScanBin(self):
        """
        int : scan bin
        """
        return self.super.scan_bin
    
    @ScanBin.setter
    def ScanBin(self, value):
        self.super.scan_bin = value

    @property
    def ChunkSize(self):
        """
        int : chunk size
        """
        return self.super.chunksize
    
    @ChunkSize.setter
    def ChunkSize(self, value):
        self.super.chunksize = value

    @property
    def DwellTime(self):
        """
        int : dwell time
        """
        return self.super.dt
    
    @DwellTime.setter
    def DwellTime(self, value):
        self.super.dt = value

    @property
    def CountImage(self):
        """
        2D numpy array : image of the total number of counts in each probe position
        """
        return np.array(self.super.Dose_image).reshape(self.super.nx,self.super.ny)


    def Run(self):
        """
        Run the 4D dataset creation
        """

        self.super.allocate_chunk() # allocate the memory for the 4D array
        self.super.init_4D_file() # initialize the hdf5 file writing
        self.super.run()

