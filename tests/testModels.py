import sys
sys.path.insert(0, "../models/model_events/build/swig")
import amici
import model_events

model = model_events.getModel()
solver = model.getSolver()


amici.readModelDataFromHDF5("./cpputest/expectedResults.h5",model.get(),"/model_events/nosensi/options")
amici.readSolverSettingsFromHDF5("./cpputest/expectedResults.h5",solver.get(),"/model_events/nosensi/options")

edata = amici.readSimulationExpData("./cpputest/expectedResults.h5","/model_events/nosensi/data",model.get())

rdata = amici.runAmiciSimulation(solver.get(),edata.get(),model.get())
print(list(rdata.x))

