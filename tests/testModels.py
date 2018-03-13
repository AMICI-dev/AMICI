import sys
sys.path.insert(0, "../build/swig/python/")
sys.path.insert(0, "../models/model_events/build")
import amici
import model_events

def trace(frame, event, arg):
    print("%s, %s:%d" % (event, frame.f_code.co_filename, frame.f_lineno))
    return trace


#sys.settrace(trace)

model = model_events.getModel()
solver = model.getSolver()

print(model.get())

amici.readModelDataFromHDF5("./cpputest/expectedResults.h5",model,"/model_events/nosensi/options")
amici.readSolverSettingsFromHDF5("./cpputest/expectedResults.h5",solver.get(),"/model_events/nosensi/options")

edata = amici.readSimulationExpData("./cpputest/expectedResults.h5","/model_events/nosensi/data",model.get())

rdata = amici.runAmiciSimulation(solver.get(),None,model.get())
print(rdata)

