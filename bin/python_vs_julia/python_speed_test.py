
import random



def make_N_random(N):
    i = 0
    while i < N:
        random.random()
        i += 1


def random_sum(N):
    i = 0
    rsum = 0
    while i < N:
        rsum += random.random() - 0.5
        i += 1
    return rsum


def fibonacchi(N):
    ia = 0
    ib = 1
    i = 0
    while i < N:
        ia, ib, i = ib, ib+ia, i+1
    return ib


if __name__ == '__main__':
    million = 1000000
    N = 100 * million
    #input("Go?")
    if False:
        print("Making %s random numbers..." % N)
        make_N_random(N)
    if False:
        print("Summing %s random numbers..." % N)
        rsum = random_sum(N)
        print("Sum: %s" % rsum)
    if True:
        N = 100
        n = input("Input number: ")
        while len(n) > 0:
            N = int(n)
            print("Calculating fibonacchi for N=%s" % N)
            fibN = fibonacchi(N)
            print(" - %s" % fibN)
            print(" - %0.04g" % fibN)
            n = input("Input number: ")
