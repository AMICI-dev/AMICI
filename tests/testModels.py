import sys
sys.path.insert(0, "../build/swig/python/")
sys.path.insert(0, "../models/model_events/build")
import amici
import model_events

def trace(frame, event, arg):
    print("%s, %s:%d" % (event, frame.f_code.co_filename, frame.f_lineno))
    return trace


sys.settrace(trace)

model = model_events.getModel()
print(model)
print(model.nx)
solver = model.getSolver()
print(solver)
rdata = amici.runAmiciSimulation(solver,None,model)
print(rdata)

