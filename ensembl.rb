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
      # Below added manually when
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
      'jgi\|Lotgi'=>['Limpet', 'Lottia gigantea']
    }
  end

  class Misc
    MISC_SPECIES_HASH = {
      /^Aqu1\./=>['Sponge', 'Amphimedon queenslandica'],
      /^AT\d/ => ['Arabidopsis','Arabidopsis thaliana']
    }
  end
end
