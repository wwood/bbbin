# A library to help make strings unique, because many bioinformatic
# programs cannot handle 
class Uniq
  def initialize()
    @names_taken ||= {}
  end
  
  def make_unique(string)
    number = 0
    to_return = string #gets changed iff it is not unique.
    while @names_taken[to_return]
      number += 1
      to_return = "#{string}-#{number}";
    end
    @names_taken[to_return] = true
    return to_return
  end
end
