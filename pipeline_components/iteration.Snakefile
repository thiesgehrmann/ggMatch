
# Config has format
# {
#   blastdbs: { genome1 : "blastDBFile" },
#   prots: { genome1 : "genome1prots.fasta" },
#   queries: "queryFile",
#   outdir: "outdir",
#   final: "finalFile"
# }

import utils as utils
import blastutils as butils


import json

###############################################################################

__OUTDIR__ = config["outdir"]

###############################################################################

rule all:
  input:
    file = config["final"]

###############################################################################

rule runBlast:
  input:
    db = lambda wildcards: config["blastdbs"][wildcards.genome],
    query = config["queries"]
  output:
    res = "%s/blastres.{genome}.tsv" % __OUTDIR__
  params:
    blastfields = butils.blastfields,
    evalue = config["evalue"]
  conda: "conda.yaml"
  shell: """
    #blastp -query {input.query} -db {input.db} -outfmt "6 {params.blastfields}" -out {output.res} -num_threads 1 -max_target_seqs 1 -evalue {params.evalue}
    diamond blastp -d {input.db} -q {input.query} -o {output.res} --max-target-seqs 1 --evalue {params.evalue} &> /dev/null
  """

###############################################################################

rule iterationFindBest:
  input:
    res   = lambda wildcards: expand("%s/blastres.{genome}.tsv" % (__OUTDIR__), genome=sorted([ g for g in config["prots"].keys() if g != "NoGenome"])),
    matchedList = config["matchedList"]
  output:
    newMatches = "%s/newMatches.tsv" % (__OUTDIR__)
  run:
    mlOrig = utils.readColumnFile(input.matchedList, "queryName iteration origin target keep", types="str int str str int")
    mlOrig = [ m for m in mlOrig if (m.keep == 1) and (int(m.iteration) > 0) ]
    mlGrouped = utils.indexListBy(mlOrig, key=lambda x: x.target.split(":")[0]).items()
    ml = { genome : set([ ":".join(g.target.split(":")[1:]) for g in group ]) for (genome, group) in mlGrouped }

    newml = []
    for (genome, resFile) in zip(sorted([ g for g in config["prots"].keys() if g != "NoGenome"]), input.res):
      bRes = butils.readBlastFile(resFile, butils.diamondfields)
        # Remove genes that have already been matched before
      bRes = [ h for h in bRes if (genome not in ml) or (h.sseqid not in ml[genome]) ]

      if len(bRes) == 0:
        continue
      #fi

      #Initially just take the best hit, if there is any
      best = butils.bestHit(bRes)
      
      newml.append( (best.qseqid, "%s:%s" % (genome,  best.sseqid)) )
    #efor

    with open(output.newMatches, "w") as ofd:
      for (origin,target) in newml:
        ofd.write("%s\t%s\n" % (origin,target))

    
###############################################################################

rule makeJSONForReciprocalBlast:
  input:
    newMatches = rules.iterationFindBest.output.newMatches
  output:
    json = "%s/recipBlast.json" % __OUTDIR__
  run:
    rbJson = { "args" : {},
               "reciprocalDir" : config["reciprocalDir"],
               "final" : "%s/reciprocalList.tsv" % __OUTDIR__,
               "evalue" : config["evalue"] }

    newMatches = utils.readColumnFile(input.newMatches, "origin target")

    for newMatch in newMatches:
      originGenome = newMatch.origin.split(":")[0]
      targetGenome = newMatch.target.split(":")[0]
      originBlastDB = config["allblastdbs"][originGenome]
      targetBlastDB = config["allblastdbs"][targetGenome]
      originQuery   = config["allprots"][originGenome]
      targetQuery   = config["allprots"][targetGenome]

      rbJson["args"]["%s:%s" % (originGenome, targetGenome)] = { "query" : originQuery, "db": targetBlastDB }
      rbJson["args"]["%s:%s" % (targetGenome, originGenome)] = { "query" : targetQuery, "db": originBlastDB }
    #efor

    with open(output.json, "w") as ofd:
      ofd.write(json.dumps(rbJson, indent=4))
    #ewith

rule runReciprocalBlast:
  input:
    json = rules.makeJSONForReciprocalBlast.output.json
  output:
    rbList = "%s/reciprocalList.tsv" % __OUTDIR__
  threads: 99
  params:
    pc_dir = config["pcDir"]
  shell: """
    {params.pc_dir}/run_reciprocal_blast.sh "{input.json}" {threads} &> /dev/null
  """


rule verifyBestFound:
  input:
    recipBlast = rules.runReciprocalBlast.output.rbList,
    newMatches = rules.iterationFindBest.output.newMatches,
    matchedList = config["matchedList"]
  output:
    fas = "%s/newOrthologs.fasta" % (__OUTDIR__)
  run:
    matches = utils.readColumnFile(input.matchedList, "query iteration origin target keep")
    newMatches = utils.readColumnFile(input.newMatches, "origin target")
    recipBlastFiles = dict([ (r.genomepair, r.file) for r in utils.readColumnFile(input.recipBlast, "genomepair file")])

    alreadyMatchedGenes = utils.indexListBy([ utils.splitgg(m.origin) for m in matches if m.keep == 1] + [ utils.splitgg(m.target) for m in matches if m.keep == 1], lambda x: x[0])
    alreadyMatchedGenes = { genome : set([g[1] for g in genes ]) for (genome, genes) in alreadyMatchedGenes.items() }

    filteredMatches = []

    for m in newMatches:
      originGenome, originGene = utils.splitgg(m.origin)
      targetGenome, targetGene = utils.splitgg(m.target)

      # Load the relevant BLAST results
      bresOT = butils.readBlastFile(recipBlastFiles.get("%s:%s" % (originGenome, targetGenome),"/dev/null"), butils.diamondfields)
      bresTO = butils.readBlastFile(recipBlastFiles.get("%s:%s" % (targetGenome, originGenome),"/dev/null"), butils.diamondfields)

      # Get the list of proteins in these genomes that have already been matched
      originAlreadyMatchedGenes = alreadyMatchedGenes.get(originGenome,set([])) - set([originGene])
      targetAlreadyMatchedGenes = alreadyMatchedGenes.get(targetGenome,set([])) - set([targetGene])

      # Remove hits to proteins that have already been matched
      bresOT = [ hit for hit in bresOT if (hit.qseqid not in originAlreadyMatchedGenes) and (hit.sseqid not in targetAlreadyMatchedGenes) ]
      bresTO = [ hit for hit in bresTO if (hit.qseqid not in targetAlreadyMatchedGenes) and (hit.sseqid not in targetAlreadyMatchedGenes) ]

      # Find the best hit per query
      bresOT = butils.bestHitPerQuery(bresOT)
      bresTO = butils.bestHitPerQuery(bresTO)

      # 
      bestForOrigin = targetGenome + ':' + ("" if originGene not in bresOT else bresOT[originGene].sseqid)
      bestForTarget = originGenome + ':' + ("" if targetGene not in bresTO else bresTO[targetGene].sseqid)

      print("########")
      print("M:  %s -> %s" % (m.origin, m.target))
      print("OT: %s -> %s" % (m.target, bestForOrigin))
      print("TO: %s -> %s" % (m.origin, bestForTarget))

      if originGenome == "__INPUTQUERY__" or ((m.origin == bestForTarget) and (m.target == bestForOrigin)):
        print(" YEEEEESSSS")
        filteredMatches.append((m.origin, m.target, True))
      else:
        print(" NOOOOOO")
        filteredMatches.append((m.origin, m.target, False))
      #fi

     # Append the matchedList
    with open(input.matchedList, "a") as ofd:
      for (origin,target, keep) in filteredMatches:
        ofd.write("%s\t%d\t%s\t%s\t%d\n" % (config["queryName"], config["iteration"], origin, target, 1 if keep else 0))
      #efor
    #ewith

      # Write the new Query file
    F = utils.loadFasta(config["queries"])
    for (origin, target, keep) in filteredMatches:
      if keep:
        targetGenome = target.split(":")[0]
        targetGene   = ":".join(target.split(":")[1:])
        prots = utils.loadFasta(config["prots"][targetGenome])
        F[target] = prots[targetGene]
      #fi
    #efor
    utils.writeFasta(F.items(), output.fas)

###############################################################################

rule final:
  input:
    file = rules.verifyBestFound.output.fas
  output:
    file = config["final"]
  shell: """
    mv {input.file} {output.file}
  """

###############################################################################
