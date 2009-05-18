# To change this template, choose Tools | Templates
# and open the template in the editor.

module Bio
  class Ensembl
    # Ensembl species taken from ensembl/modules/Bio/EnsEMBL/Registry.pm
    # Ensembl version 52 (December 2009)
    # transformed manually to hash form
    ENSEMBL_SPECIES_HASH = {'ENSRNO'=>['Rat','Rattus norvegicus'],
      'ENSMUS'=>['Mouse','Mus musculus'],
      'ENSGAL'=>['Chicken','Gallus gallus'],
      'ENSBTA'=>['Cow','Bos taurus'],
      'ENSDAR'=>['Zebrafish','Danio rerio'],
      'ENSCAF'=>['Dog','Canis familiaris'],
      'ENSPTR'=>['Chimpanzee','Pan troglodytes'],
      'ENSCPO'=>['Guinea pig','Cavia porcellus'],
      'ENSCIN'=>['C. intestinalis','Ciona intestinalis'],
      'ENSCSAV'=>['C. savignyi','Ciona savignyi'],
      'ENSDNO'=>['Armadillo','Dasypus novemcinctus'],
      'ENSETE'=>['Lesser hedgehog tenrec','Echinops telfairi'],
      'ENSEEU'=>['Hedgehog','Erinaceus europaeus'],
      'ENSFCA'=>['Cat','Felis catus'],
      'ENSGAC'=>['Stickleback','Gasterosteus aculeatus'],
      'ENSLAF'=>['Elephant','Loxodonta africana'],
      'ENSMMU'=>['Macaque','Macaca mulatta'],
      'ENSMOD'=>['Opossum','Monodelphis domestica'],
      'ENSMLU'=>['Microbat','Myotis lucifugus'],
      'ENSOAN'=>['Platypus','Ornithorhynchus anatinus'],
      'ENSOCU'=>['Rabbit','Oryctolagus cuniculus'],
      'ENSORL'=>['Medaka','Oryzias latipes'],
      'ENSSAR'=>['Shrew','Sorex araneus'],
      'ENSSTO'=>['Squirrel','Spermophilus tridecemlineatus'],
      'ENSTBE'=>['Tree shrew','Tupaia belangeri'],
      'SINFRU'=>['Fugu','Takifugu rubripes'],
      'ENSXET'=>['Frog','Xenopus tropicalis'],
      # Below added manually when unknown genes are found.
      'ENSPCA'=>['Hyrax','Procavia capensis'],
      'ENSOPR'=>['Pika','Ochotona princeps'],
      'ENSDOR'=>['Kangaroo rat','Dipodomys ordii'],
      'ENSECA'=>['Horse','Equus caballus'],
      'ENSPVA'=>['Megabat','Pteropus vampyrus'],
      'ENSTTR'=>['Dolphin','Tursiops truncatus'],
      'ENSMIC'=>['Mouse lemur','Microcebus murinus'],
      'ENSOGA'=>['Bushbaby','Otolemur garnettii'],
      'ENSTRU'=>['Fugu','Takifugu rubripes'],
      'ENSTNI'=>['Tetraodon','Tetraodon nigroviridis'],
      'ENSTSY'=>['Tarsier','Tarsius syrichta'],
      'ENSVPA'=>['Alpaca','Vicugna pacos'],
      'ENSPPY'=>['Orangutan','Pongo pygmaeus'],
      'ENSGGO'=>['Gorilla','gorilla gorilla'],
      'ENSACA'=>['Anole lizard','Anolis carolinensis'],
      'ENSTGU'=>['Zebra finch','Taeniopygia guttata'],
      'ENS'=>['Human','Homo sapiens']
    }
    ENSEMBL_OTHER_HASH = {
      /^Y[A-Z][A-Z]\d/ => ['Yeast','Saccharomyces cerevisiae'],
      /^T\d/ => ['Worm','Caenorhabditis elegans'],
      /^F\d/ => ['Worm','Caenorhabditis elegans'],
      /^Y\d/ => ['Worm','Caenorhabditis elegans'],
      /^FBpp\d/ => ['Fly','Drosophila melanogaster'],
      /^AAEL\d/ => ['Aedes','Aedes aegypti'],
      /^AGAP\d/ => ['Anopholes','Anopholes gambiae']
    }
  end

  class JGI
    JGI_SPECIES_HASH = {
      'jgi\|Brafl'=>['Lancelet', 'Branchiostoma floridae'],
      'jgi\|Helro'=>['Leech', 'Helobdella robusta'],
      'jgi\|Triad'=>['Trichoplax', 'Trichoplax adhaerens'],
      'jgi\|Monbr'=>['Choanoflagellate', 'Monosiga brevicollis'],
      'jgi\|Nemve'=>['Sea anemone', 'Nematostella vectensis'],
      'jgi\|Lotgi'=>['Limpet', 'Lottia gigantea'],
      'jgi\|Thaps'=>['Diatom', 'Thalassiosira pseudonana'],
      'jgi\|Chlre'=>['Green Alga', 'Chlamydomonas reinhardtii']
    }
  end
  
  class EuPathDB
    EUPATHDB_SPECIES_HASH = {
      /gb\|GL.* organism=Giardia_lamblia/ => ['Giardia', 'Giardia lamblia']
    }
  end

  class Misc
    MISC_SPECIES_HASH = {
      /^Aqu1\./=>['Sponge', 'Amphimedon queenslandica'],
      /^AT\d/ => ['Arabidopsis','Arabidopsis thaliana'],
      /^gnl\|CMER/ => ['Red Alga','Cyanidioschyzon merolae'],
      /\|\|\|\|Entamoeba histolytica\|/ => ['Entamoeba','Entamoeba histolytica']
    }
  end

  class Orthomcl
    ORTHOMCL_SPECIES_HASH = {
/^aae\|/ => ['Aquifex aeolicus VF5'],
/^tma\|/ => ['Thermotoga maritima MSB8'],
/^det\|/ => ['Dehalococcoides ethenogenes 195'],
/^dra\|/ => ['Deinococcus radiodurans R1'],
/^tpa\|/ => ['Treponema pallidum subsp. pallidum str. Nichols'],
/^cte\|/ => ['Chlorobium tepidum TLS'],
/^rba\|/ => ['Rhodopirellula baltica SH 1'],
/^cpn\|/ => ['Chlamydophila pneumoniae CWL029'],
/^syn\|/ => ['Synechococcus sp. WH 8102'],
/^mtu\|/ => ['Mycobacterium tuberculosis H37Rv'],
/^ban\|/ => ['Bacillus anthracis str. Ames Ancestor'],
/^wsu\|/ => ['Wolinella succinogenes DSM 1740'],
/^gsu\|/ => ['Geobacter sulfurreducens PCA'],
/^atu\|/ => ['Agrobacterium tumefaciens str. C58'],
/^rso\|/ => ['Ralstonia solanacearum GMI1000'],
/^eco\|/ => ['Escherichia coli W3110'],
/^ftu\|/ => ['Francisella tularensis subsp. tularensis SCHU S4'],
/^ype\|/ => ['Yersinia pestis CO92'],
/^sfl\|/ => ['Shigella flexneri 2a str. 301'],
/^sty\|/ => ['Salmonella enterica subsp. enterica serovar Typhi str. CT18'],
/^sau\|/ => ['Staphylococcus aureus subsp. aureus Mu50'],
/^vch\|/ => ['Vibrio cholerae O1 biovar eltor str. N16961'],
/^lmo\|/ => ['Listeria monocytogenes EGD-e'],
/^cje\|/ => ['Campylobacter jejuni subsp. jejuni NCTC 11168'],
/^spn\|/ => ['Streptococcus pneumoniae TIGR4'],
/^cpe\|/ => ['Clostridium perfringens str. 13'],
/^bur\|/ => ['Burkholderia mallei ATCC 23344'],
/^rty\|/ => ['Rickettsia typhi str. Wilmington'],
/^bsu\|/ => ['Brucella suis 1330'],
/^cbu\|/ => ['Coxiella burnetii RSA 493'],
/^hal\|/ => ['Halobacterium sp. NRC-1'],
/^mja\|/ => ['Methanocaldococcus jannaschii DSM 2661'],
/^sso\|/ => ['Sulfolobus solfataricus P2'],
/^neq\|/ => ['Nanoarchaeum equitans Kin4-M'],
/^gla\|/ => ['Giardia lamblia ATCC 50803'],
/^ehi\|/ => ['Entamoeba histolytica HM-1:IMSS'],
/^ddi\|/ => ['Dictyostelium discoideum AX4'],
/^pfa\|/ => ['Plasmodium falciparum 3D7'],
/^pyo\|/ => ['Plasmodium yoelii yoelii str. 17XNL'],
/^pvi\|/ => ['Plasmodium vivax SaI-1'],
/^pkn\|/ => ['Plasmodium knowlesi strain H'],
/^pbe\|/ => ['Plasmodium berghei'],
/^pch\|/ => ['Plasmodium chabaudi'],
/^cpa\|/ => ['Cryptosporidium parvum'],
/^cho\|/ => ['Cryptosporidium hominis'],
/^tgo\|/ => ['Toxoplasma gondii'],
/^the\|/ => ['Theileria parva'],
/^tan\|/ => ['Theileria annulata'],
/^tth\|/ => ['Tetrahymena thermophila'],
/^tcr\|/ => ['Trypanosoma cruzi'],
/^tbr\|/ => ['Trypanosoma brucei'],
/^lma\|/ => ['Leishmania major'],
/^ecu\|/ => ['Encephalitozoon cuniculi'],
/^sce\|/ => ['Saccharomyces cerevisiae'],
/^spo\|/ => ['Schizosaccharomyces pombe'],
/^cne\|/ => ['Filobasidiella neoformans'],
/^pha\|/ => ['Phanerochaete chrysosporium'],
/^ago\|/ => ['Eremothecium gossypii'],
/^ncr\|/ => ['Neurospora crassa'],
/^aor\|/ => ['Aspergillus oryzae'],
/^yli\|/ => ['Yarrowia lipolytica'],
/^kla\|/ => ['Kluyveromyces lactis'],
/^dha\|/ => ['Debaryomyces hansenii'],
/^cgl\|/ => ['Candida glabrata'],
/^gth\|/ => ['Guillardia theta'],
/^cre\|/ => ['Chlamydomonas reinhardtii'],
/^ota\|/ => ['Ostreococcus tauri'],
/^cme\|/ => ['Cyanidioschyzon merolae'],
/^tps\|/ => ['Thalassiosira pseudonana'],
/^ath\|/ => ['Arabidopsis thaliana'],
/^osa\|/ => ['Oryza sativa'],
/^cel\|/ => ['Caenorhabditis elegans'],
/^cbr\|/ => ['Caenorhabditis briggsae'],
/^bma\|/ => ['Brugia malayi'],
/^dme\|/ => ['Drosophila melanogaster'],
/^aga\|/ => ['Anopheles gambiae'],
/^sma\|/ => ['Schistosoma mansoni'],
/^aed\|/ => ['Aedes aegypti'],
/^ame\|/ => ['Apis mellifera'],
/^cin\|/ => ['Ciona intestinalis'],
/^fru\|/ => ['Takifugu rubripes'],
/^tni\|/ => ['Tetraodon nigroviridis'],
/^dre\|/ => ['Danio rerio'],
/^gga\|/ => ['Gallus gallus'],
/^mmu\|/ => ['Mus musculus'],
/^rno\|/ => ['Rattus norvegicus'],
/^hsa\|/ => ['Homo sapiens']
}
  end
end
