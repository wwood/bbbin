class Hydrophobicity
  HYDROPHOBICITY = {
      'A' => 1.8,
      'R' => -4.5,
      'N' => -3.5,
      'D' => -3.5,
      'C' => 2.5,
      'Q' => -3.5,
      'E' => -3.5,
      'G' => -0.4,
      'H' => -3.2,
      'I' => 4.5,
      'L' => 3.8,
      'K' => -3.9,
      'M' => 1.9,
      'F' => 2.8,
      'P' => -1.6,
      'S' => -0.8,
      'T' => -0.7,
      'W' => -0.9,
      'Y' => -1.3,
      'V' => 4.2
  }
  
  # Taken from http://en.wikipedia.org/wiki/Proteinogenic_amino_acid
  # July 8, 2010
  PKA = {
  'D' => 3.9,
  'E' => 4.07,
  'H' => 6.04,
  'K' => 10.54,
  'R' => 12.48,
  'U' => 5.73,
  'Y' => 10.46,
  }
  
  # Given an amino acid sequence, return an array of
  # hydrophobicity values. They may contain nil
  # values if an X was encountered, for instance
  def hydrophobicity_profile(amino_acid_sequence)
    profile(amino_acid_sequence, HYDROPHOBICITY)
  end
  
  # Return a profile (e.g. hydrophobicity) 
  # for each amino acid in the sequence, where the
  # profile hash is a hash of single letter amino
  # acid names to values of the profile. The profile_hash
  # need not account for each amino acid in the sequence.
  # In that case a nil is returned as part of the profile
  def profile(amino_acid_sequence, profile_hash)
    acmi = []
    amino_acid_sequence.each_char do |amino|
      if profile_hash[amino]
        acmi.push profile_hash[amino]
      else
        acmi.push nil # not recorded in the scale
      end
    end
    acmi
  end
end