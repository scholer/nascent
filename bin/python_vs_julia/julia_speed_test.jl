

function make_N_random(N::Integer)
  i = 0
  while i < N
    rand()
    i += 1
  end
end

function random_sum(N::Integer)
  i = 0
  sum = 0
  while i < N
    sum += rand() - 0.5
    i += 1
  end
  return sum
end

function fibonacchi(N::Integer)
  ia = Int128(0)
  ib = Int128(1)
  i = Int128(0)
  while i < N
    #println("-\n $ia, $ib, $i")
    ia, ib, i = ib, ib+ia, i+1
    #println(" $ia, $ib, $i")
    # 1, 1
    # 1, 2
    # 2, 3
    # 3, 5
    # 5, 8
  end
  return ib
end


function main()
  println("Starting...")
  million = 1000000
  N = 100 * million
  #input("Go? ")
  #print("Making $N random numbers...")
  #make_N_random(N) # This is simply skipped - optimized away by the compiler.
  println("Done!")
  print("Summing $N random numbers...")
  @time rsum = random_sum(N)
  println("Sum: $rsum")
  #N = 100
  #n = input("Input number: ")
  n = "92"
  while length(n) > 0
    N = parse(Int, n)
    println("Calculating fibonacchi for N=$N")
    fibN = fibonacchi(N)
    println(" - $fibN - ($(1.0*fibN))")
    println(" - type: $(typeof(fibN))")
    n = input("\nInput number: ")
  end
end

function input(prompt::String="")
  ## Note: Pressing ctrl+C during input can cause the program to halt
  print(prompt)
  chomp(readline())
end



println("Entering main()...")
main()

#=
When it comes to random number generation, pypy and julia are about on par. CPython is maybe 10 times slower.


=#
