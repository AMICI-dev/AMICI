#!/usr/bin/env python3

import sys
import numpy as np

from example import AmiciExample

class ExampleNeuron(AmiciExample):

    def __init__(self):
        AmiciExample.__init__( self )

        self.numZ = 1
        self.numX = 2
        self.numP = 4
        self.numK = 2

        self.modelOptions['theta'] = np.log10([0.02, 0.3, 65, 0.9])
        self.modelOptions['kappa'] = [-60,10]
        self.modelOptions['ts'] = [0.0, 0.9991997694102649, 1.9983995388205298, 2.9975993082307943, 3.9967990776410596, 4.995998847051324, 5.995198616461589, 6.994398385871854, 7.993598155282119, 8.992797924692384, 9.991997694102649, 10.991197463512913, 11.990397232923177, 12.989597002333443, 13.988796771743708, 14.987996541153972, 15.987196310564238, 16.986396079974504, 17.98559584938477, 18.984795618795033, 19.983995388205297, 20.98319515761556, 21.982394927025826, 22.981594696436094, 23.980794465846355, 24.979994235256623, 25.979194004666887, 26.978393774077155, 27.977593543487416, 28.976793312897684, 29.975993082307944, 30.975192851718212, 31.974392621128477, 32.97359239053874, 33.97279215994901, 34.97199192935927, 35.97119169876954, 36.9703914681798, 37.969591237590066, 38.96879100700033, 39.967990776410595, 40.967190545820856, 41.96639031523112, 42.96559008464139, 43.96478985405165, 44.96398962346192, 45.96318939287219, 46.962389162282456, 47.96158893169271, 48.96078870110298, 49.959988470513245, 50.959188239923506, 51.958388009333774, 52.95758777874404, 53.95678754815431, 54.95598731756456, 55.95518708697483, 56.9543868563851, 57.95358662579537, 58.95278639520563, 59.95198616461589, 60.95118593402616, 61.950385703436424, 62.949585472846685, 63.94878524225695, 64.94798501166721, 65.94718478107748, 66.94638455048775, 67.94558431989802, 68.94478408930829, 69.94398385871854, 70.9431836281288, 71.94238339753907, 72.94158316694934, 73.9407829363596, 74.93998270576986, 75.93918247518013, 76.9383822445904, 77.93758201400065, 78.93678178341092, 79.93598155282119, 80.93518132223146, 81.93438109164171, 82.93358086105198, 83.93278063046225, 84.93198039987251, 85.93118016928278, 86.93037993869305, 87.9295797081033, 88.92877947751359, 89.92797924692384, 90.9271790163341, 91.92637878574438, 92.92557855515463, 93.92477832456491, 94.92397809397517, 95.92317786338542, 96.9223776327957, 97.92157740220595, 98.92077717161622]
        self.modelOptions['pscale'] = 2

        self.solverOptions['atol'] = 1e-16
        self.solverOptions['maxsteps'] = 1e5
        self.solverOptions['nmaxevent'] = 22
        self.solverOptions['rtol'] = 1e-12
        self.solverOptions['sens_ind'] = []
        self.solverOptions['sensi'] = 0
        self.solverOptions['sensi_meth'] = 1

        self.data['Y'] = np.full((len(self.modelOptions['ts']), 1), np.nan)
        self.data['Sigma_Y'] = np.full((len(self.modelOptions['ts']), 1), np.nan)

        self.data['Z'] = np.transpose([[2.4420740245701733, 5.921007525639647, np.nan, 10.366527794075543, 12.83694308382395, 14.624269559247253, 18.722446363647578, 22.20739095602005, 28.602369747827655, 31.442843729542822, 34.01927181474919, 41.26726577405225, 44.275254172160395, 51.56486254598814, 57.100273114298204, 61.961654997481084, 69.03838073191332, 74.3546047146856, 81.21960802401809, 87.2873927650102, 93.34894804384085, 98.57346300859241]])
        self.data['Sigma_Z'] = 0.5 * np.ones(self.data['Z'].shape)

        self.data['condition'] = self.modelOptions['kappa']
        self.data['t'] = self.modelOptions['ts']

    """   
        TODO: should generate from simulation: 
        rng(0);
        D.Z = sol.z(sol.z<t(end));
        t = linspace(0,D.Z(end)-0.1,100);
        D.Z = D.Z + 0.5*randn(size(D.Z));
        D.Z(3) = NaN;
        D.Sigma_Z = 0.5*ones(size(D.Z));
        D.Z = D.Z + D.Sigma_Z.*randn(size(D.Z));       
    """


def writeNoSensi(filename):
    ex = ExampleNeuron()

    ex.writeToFile(filename, '/model_neuron/nosensi/')


def writeSensiForward(filename):
    ex = ExampleNeuron()

    ex.solverOptions['sens_ind'] = np.arange(0, ex.numP)
    ex.solverOptions['sensi'] = 1

    ex.writeToFile(filename, '/model_neuron/sensiforward/')


def writeSensi2Forward(filename):
    ex = ExampleNeuron()

    ex.solverOptions['sens_ind'] = np.arange(0, ex.numP)
    ex.solverOptions['sensi'] = 2
    ex.solverOptions['sensi_meth'] = 1

    ex.writeToFile(filename, '/model_neuron/sensi2forward/')

def main():
    if len(sys.argv) < 2:
        print("Error: Must provide output file as first and only argument.")
        sys.exit(1)
    filename = sys.argv[1]

    writeNoSensi(filename)
    writeSensiForward(filename)
    writeSensi2Forward(filename)

if __name__ == "__main__":
    main()
