# FAQ

__Q__: My model fails to build.

__A__: Remove the corresponding model directory located in AMICI/models/*yourmodelname* and compile again.

---

__Q__: It still does not compile.

__A__: Remove the directory AMICI/models/`mexext` and compile again.

---

__Q__: It still does not compile.

__A__: Make an [issue](https://github.com/ICB-DCM/AMICI/issues) and we will have a look.

---

__Q__: My Python-generated model does not compile from MATLAB.

__A__: Try building any of the available examples before. If this succeeds, 
retry building the original model. Some dependencies might not be built 
correctly when using only the `compileMexFile.m` script. 

---

__Q__: I get an out of memory error while compiling my model on a Windows machine.

__A__: This may be due to an old compiler version. See [issue #161](https://github.com/ICB-DCM/AMICI/issues/161) for instructions on how to install a new compiler.

---

__Q__: How are events interpreted in a DAE context?

__A__: Currently we only support impulse free events. Also sensitivities have never been tested. Proceed with care and create an [issue](https://github.com/ICB-DCM/AMICI/issues) if any problems arise!

---

__Q__: The simulation/sensitivities I get are incorrect.

__A__: There are some known issues, especially with adjoint sensitivities, events and DAEs. If your particular problem is not featured in the [issues](https://github.com/ICB-DCM/AMICI/issues) list, please add it!

