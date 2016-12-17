# An oscillatory neural network model of motor dynamics during continuous periodic movement

### authors:

Wisam Reid<br />
Center for Computer Research in Music and Acoustics<br />
Stanford University<br />
wisam@ccrma.stanford.edu<br />

Iran Roman<br />
Stanford Neuroscience Graduate Program<br />
Center for Computer Research in Music and Acoustics<br />
Stanford University<br />
iran@stanford.edu<br />

### abstract: 

Brain oscillations are of functional relevance to the healthy
functioning of the brain. In the motor system, desynchronization of
oscillatory activity in the beta band ( 20Hz) is observed through
periodic somatosensory and/or auditory stimulation. This
desynchronization occurs at a rate equal to the period of stimulation,
and resynchronization of oscillatory activity anticipates the next
stimulus. The computations underlying motor function can be explained by
networks of neural oscillators carrying out non-linear transformations
of stimuli, reflected by desynchronization of oscillatory activity in
the motor system. We developed a neural model of nonlinear oscillators
that transforms stimuli into characteristic oscillatory activity of the
beta band. Our model is built using a canonical model of nonlinear
oscillators, captures even-related desynchronization in the beta band,
and is able to anticipate the next period of stimulation through
resynchronization. Additionally, our model provides a dynamical-systems
perspective of computations and stimuli transformations in the motor system.

-----  

## The Code
### Fixed Point Analysis and Amplitude Vector Plane Plots 
#### [Folder] analysis/

      script in this folder utilize functions in lib/

##### autonomousOscAnalysis.m

      Autonomous oscillator parameter analysis:      
      Plots the amplitude vector feild for ranging oscillator parameters
      
      Current Total Runtime: ~70 seconds for 9 images in the 
      12-16-2016 research report

### Models and Demos
#### [Folder] models/

      Models and prototypes 
      
#### [Folder] demos/

      1. Demo for McClelland Lab Meeting
      2. Bursting Oscillator 3D animation 
      
#### [Folder] images/

      Matlab Images 
      
### Functions

#### [Folder] lib/

      type "help name_of_file" in the Matlab terminal for more information 
            1. getFP **Passes 66 Unit Test** (tests were made from figures in Kim & Large 2015)
            2. getSS **Passes the same 66 Unit Test** in tests/Fixed_Points
            3. getOscParamRegime **Passes 10 Unit Test** in tests/Param_Regime
            4. plotAmplitudeVectorFeild
            5. plotFP

### Unit Tests
#### [Folder] tests/

            To run tests use  
            1. runFixedPointUnitTests
            2. runOscParamRegimeUnitTests
            
-----
