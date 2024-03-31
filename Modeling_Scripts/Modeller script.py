from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()  # Enable verbose output
env = environ() # Initialize the Modeller environment

# Directory where the input data is stored and where the output will go
env.io.atom_files_directory = ['./', 'Models']
aln = alignment(env)
aln.read(file= 'alignment.pir' ,alignment_format= 'PIR') 
# Build the homology models
a = loopmodel(env,
            alnfile='alignment.pir',  # PIR alignment filename
            knowns='1C8Q_clean',      # code of the template
            sequence='AMY_PSEHA',     # code of the target
            assess_methods=(assess.DOPE,     
                            assess.GA341))    
a.starting_model = 1
a.ending_model = 10  # Number of models to build
a.make() #  Build the models