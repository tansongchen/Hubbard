from scipy.special import comb

s = sum([comb(8,i)**2 for i in range(9)])

print(s)