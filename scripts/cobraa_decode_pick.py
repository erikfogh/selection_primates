import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='comma-separated list of paths', nargs="+")
args = parser.parse_args()

def get_LL_from_file(filename):
    with open(filename) as f:
        lines = f.readlines()
    return float([i for i in lines if 'likelihood' in i and 'final' in i][0].split(' ')[-1])


best_param_file = None
best_val = None
# Go through the files and pick the 
for p in (args.i):
    LL = get_LL_from_file(p)
    if best_val == None:
        best_val = LL
        best_param_file = p
    if best_val < LL:
        best_val = LL
        best_param_file = p
    

with open(best_param_file) as f:
    finallines = f.readlines()
ztheta = float([i for i in finallines if 'theta' in i ][0].split(' ')[-1])
zrho = float([i for i in finallines if 'rho' in i ][0].split(' ')[-1])
zgamma = float([i for i in finallines if 'gamma' in i ][0].split(' ')[-1])
file_name = best_param_file.split("/")[-1]
zte = file_name.split("te")[1].split("_")[0]
zts = file_name.split("ts")[1].split("_")[0]
    
final_params = np.loadtxt(best_param_file)
lambdaA_parameters = ",".join([str(x) for x in final_params[:,2]*ztheta/4])
lambdaB_parameters = ",".join([str(x) for x in final_params[:,3]*ztheta/4])

out_print = """
-theta {ztheta} -rho {zrho} -gamma_fg {zgamma} \
-lambda_A_fg {lambdaA_parameters} -lambda_B_fg {lambdaB_parameters} \
-ts {zts} -te {zte}
""".format(ztheta=ztheta, zrho=zrho, zgamma=zgamma,
lambdaA_parameters=lambdaA_parameters, lambdaB_parameters=lambdaB_parameters,
zts=zts, zte=zte)

print(out_print)