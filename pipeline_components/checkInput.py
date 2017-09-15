#!/usr/bin/env python

import json
import sys
import inspect, os

import utils as utils

__INSTALL_DIR__ = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

if len(sys.argv) < 3:
  print("Error, incomplete arguments to checkInput.py")
  sys.exit(1)
#fi

configFile = sys.argv[1];
action     = sys.argv[2];


config = {}
if not(os.path.isfile(sys.argv[1])):
  errors.append("ConfigFile '%s'doesn't exist!"% sys.argv[1])
  config = { "genomes": {}, "queries" : {}, "outdir" : "/dev/null" }
else:
  config = json.load(open(sys.argv[1],"r"))
#fi
dconfig = json.load(open("%s/defaults.json" % __INSTALL_DIR__,"r"))

errors = []
warnings = []

if "genomes" not in config:
  errors.append("No genomes defined!")
else:
  for genome in config["genomes"]:
    if "prots" not in config["genomes"][genome]:
      errors.append("No protein file specified for genome '%s'." % genome)
    else:
      if not(os.path.isfile(config["genomes"][genome]["prots"])):
        errors.append("Protein file specified for genome '%s' does not exist." % genome)
      #fi
    #fi
  #efor
#fi

loadedGenomes = {}

if "queries" not in config:
  errors.append("No queries defined!")
else:
  for query in config["queries"]:
    if not(os.path.isfile(config["queries"][query])):
      errors.append("Query file specified for query '%s' des not exist." % genome)
    else:
      F = utils.loadFasta(config["queries"][query])
      for seed in F:
        genome, gene = utils.splitgg(seed)
        if genome is None:
          continue
        #fi
        if genome not in config["genomes"]:
          errors.append("Genome '%s' referenced in seed sequence '%s' for query '%s' is not present in genomes list" % (genome, seed, query))
        else:
          if (genome not in loadedGenomes) and ("prots" in config["genomes"][genome]) and os.path.isfile(config["genomes"][genome]["prots"]):
            loadedGenomes[genome] = utils.loadFasta(config["genomes"][genome]["prots"])
          #fi
          relGenome = loadedGenomes[genome]
          if gene not in relGenome:
            errors.append("Gene '%s' referenced in seed sequence '%s' in query '%s' is not present in genome '%s'." % (gene, seed, query, genome))
          #fi
        #else
      #efor
  #efor
#fi

if "outdir" not in config:
  warnings.append("Outdir is not specified. Defaulting to '%s'." % dconfig["outdir"])
#fi

if action == "phylogeny":
  if "phylogeny_reference" not in config:
    errors.append("Field \"phylogeny_reference\" is not defined in the configuration file.")
  else:
    if config["phylogeny_reference"] not in config["genomes"]:
      errors.append("The phylogeny reference genome \"%s\" is not present in the list of genomes." % config["phylogeny_reference"])
    #fi
  #fi

  if "phylogeny_outgroups" not in config:
    warnings.append("No phylogenetic outgroup is defined in 'phylogeny_outgroups'. Will perform midpoint re-rooting.")
  else:
    for outgroup in config['phylogeny_outgroups'].split(' '):
      if outgroup not in config["genomes"]:
        error.append("Outgroup genome '%s' specified in 'phylogeny_outgroups' is not present in the list of genomes." % outgroup)
      #fi
    #efor
  #fi
#fi

if action == "validate":
  if utils.which("interproscan.sh") is None:
    errors.append("`interproscan.sh` is not in $PATH. (PATH=%s)" % os.environ["PATH"])
  #fi
#fi



for error in errors:
  print("ERROR: %s" % error)
#efor

for warning in warnings:
  print("WARNING: %s" % warning)
#efor

if len(errors) > 0:
  sys.exit(1)
#fi

sys.exit(0)
