# Package validation

Validation program include public datasets and generated datasets. All public datasets include in test/testdata.jl.

Generated datasets made with *randrbeds* function with fixed seed (could be reproduced at any time) and include following disignes:

### 24 subjects, balanced

  * 1  TRTR/RTRT
  * 2  TRRT/RTTR
  * 3  TTRR/RRTT
  * 4  TRTR/RTRT/TRRT/RTTR
  * 5  TRRT/RTTR/TTRR/RRTT
  * 6  TRTR/RTRT/TTRR/RRTT
  * 7  TRT/RTR
  * 8  TRR/RTT
  * 9  TR/RT/TT/RR
  * 10 TRR/RTR/RRT
  * 11 TRR/RTR*

  SAS procedures for generated datasets can be found in [*validation/sas/*](https://github.com/PharmCat/ReplicateBE.jl/tree/master/validation/sas) folder.

  All validation datasets include in package test procedure.

  *SAS WARNING: Stopped because of infinite likelihood.
