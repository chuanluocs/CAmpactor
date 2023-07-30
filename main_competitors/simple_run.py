import sys, os, argparse, random
from formatencoding import FormatEncoder
from formattranslater import format_translate

fullname_dict = {"A": "AutoCCAG", "F": "FastCA", "T": "TCA"}

subdir_dict = fullname_dict

binary_path_dict = {"A": "MetaTCA/TCA", "F": "FastCA", "T": "TCA"}

# commands for running AutoCCAG, FastCA and TCA respectively
solver_cmds = {"A": "python3 run_autoccag.py {model} {constraints} {seed} {cutoff} {initca}", 
    "F": "./FastCA {model} {constraints} {cutoff} {seed} {initca}", 
    "T": "./TCA {model} {constraints} {cutoff} {seed} {initca}"
}

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="A unified interface for running PCAG solvers")
    parser.add_argument("-s", "--solver", help="choose the target PCAG solver (A for AutoCCAG, F for FastCA, T for TCA)", 
        choices=["A", "F", "T"], required=True)
    parser.add_argument("-i", "--input_CNF_file", help="path of input CNF", required=True)
    parser.add_argument("-o", "--output_PCA_path", help="path of output PCA", required=True)
    parser.add_argument("--seed", help="random seed", type=int, default=1)
    parser.add_argument("--cutoff", help="cutoff time (in seconds)", type=float, default=3600.0)
    parser.add_argument("--init_PCA_path", metavar="PATH", help="path of initial PCA (if not given, SamplingCA will be called to generate one)")

    args = parser.parse_args()

    if not os.path.exists(args.input_CNF_file):
        print("c Error: input CNF file not found. ")
        sys.exit(-1)
    
    if args.init_PCA_path:
        if not os.path.exists(args.init_PCA_path):
            print("c Error: initial PCA not found. ")
            sys.exit(-1)
    else:
        if not os.path.exists("../SamplingCA/SamplingCA"):
            print("c Error: SamplingCA binary not found. ")
            sys.exit(-1)

    if not os.path.exists(os.path.join(subdir_dict[args.solver], binary_path_dict[args.solver])):
        print("c Error: please make sure %s is built correctly. " % fullname_dict[args.solver])
        sys.exit(-1)
    
    random.seed(os.getpid())
    randstr = str(os.getpid()) + "_" + str(random.randint(1, sys.maxsize))

    real_initca = ""
    if not args.init_PCA_path:
        print("c running SamplingCA ...")
        real_initca = "%s_tmp_initial_CA.out" % randstr
        real_initca = os.path.abspath(real_initca)
        # command for running SamplingCA
        real_input_path = os.path.abspath(args.input_CNF_file)
        os.chdir("../SamplingCA")
        cmd = "./SamplingCA -seed %d -input_cnf_path %s -output_testcase_path %s" % (args.seed, real_input_path, real_initca)
        os.system(cmd)
        os.chdir("../main_competitors")
    else:
        real_initca = os.path.abspath(args.init_PCA_path)

    fe = FormatEncoder(args.input_CNF_file, 2, randstr)
    fe.encoding()

    print("c running %s ..." % (fullname_dict[args.solver]))
    os.chdir(subdir_dict[args.solver])
    cmd = solver_cmds[args.solver].format(model="../" + fe.model_file_path, constraints="../" + fe.constraints_file_path, 
        seed=args.seed, cutoff=int(args.cutoff), initca=real_initca)
    os.system("timeout %fs %s" % (args.cutoff, cmd))

    print("c done. ")
    cur_output = "___%s_temp_CA_file" % fullname_dict[args.solver]
    cur_output = os.path.join(subdir_dict[args.solver], cur_output)
    os.chdir("..")

    if os.path.exists(cur_output):
        format_translate(cur_output, args.output_PCA_path)
        os.system("rm " + cur_output)
    else:
        print("c the initial PCA is not optimized at all. ")
        os.system("cp %s %s" % (real_initca, args.output_PCA_path))

    fe.clear()
    if not args.init_PCA_path:
        os.system("rm " + real_initca)
