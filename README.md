
# EvenTem
EvenTem started as a fork of [riCOM](https://github.com/ThFriedrich/riCOM_cpp), but has grown to a complete tool for fast 4D STEM reconstruction, with a focus on **live update capabilities** and optimization for **event-driven pipelines**. EvenTem is fully written in C++, but has a python API through [pybind11](https://github.com/pybind/pybind11), so that it can be simply installed and used like any python module.

API documentation can be found [here]().

Supported detectors are
- [Amsterdam Scientific Instruments CheeTah T3](https://www.amscins.com/product/cheetah-series/) (.tpx3)
- [Advacam AdvaPIX TPX3](https://advacam.com/camera/advapix-tpx3/) (.t3p)
- [Quantum detectors MerlinEM](https://quantumdetectors.com/products/merlinem/) (.mib)

as well as data in the following formats
- Frame-based [numpy](https://numpy.org/doc/stable/index.html) files (dtype np.uint8 & np.uint16, shape [RX,RY,K,K] with K either 64,128,256,512)
- Event-based .electron (custom dtype [RX,RY,KX,XY,scan number], all uint16)
  

## Installation
Pre-compiled binaries should be available for Windows, Linux and Apple silicon supporting python 3.9 to 3.12. For these, installation requires to clone the repository using
```
git clone https://github.com/EMAT-Jo/EvenTem.git
```
and then to navigate to the eventSTEM directory and use
```
pip install .
```
Testing the installation can be done as
```python
import EvenTem
```
## GUI Use
see EvenTem GUI repo
## Scripting Use
An example reconstruction requires as little as
```python
from EvenTem import vSTEM
vs = vSTEM(nx=N,ny=N,repetitions=repetitions,filename=filename)
vs.DetectorSize = DetectorSize
vs.InnerRadia = [vs.DetectorSize//2]
vs.OuterRadia = [vs.DetectorSize]
vs.Offsets = [(vs.DetectorSize//2,vs.DetectorSize//2)]
vs.Run()
vs.PlotImage()

#results are numpy accessible
result = vs.image #np.array
```
For complete example notebooks, check the examples directory.
