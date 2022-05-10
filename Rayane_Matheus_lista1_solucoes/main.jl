include("optimizerEngine.jl")

# LETRA A #
oper = "Max"
nvar = 2
nrest = 2
indices = [
    2 1 6;
    1 3 9;
    3 1 0; #F.O
]
sinal = [
    "≤",
    "≤",
]

a = [oper, nvar, nrest, indices, sinal]

# LETRA B #
oper = "Max"
nvar = 2
nrest = 2
indices = [
    1 1 4;
    1 -1 5;
    1 1 0; #F.O
]
sinal = [
    "≤",
    "≤",
]

b = [oper, nvar, nrest, indices, sinal]

# LETRA C #
oper = "Max"
nvar = 2
nrest = 2
indices = [
    8 2 16;
    5 2 12;
    4 1 0; #F.O
]
sinal = [
    "≤",
    "≤",
]

c = [oper, nvar, nrest, indices, sinal]

# LETRA D #
oper = "Max"
nvar = 2
nrest = 2
indices = [
    1 -1 4;
    1 2 4;
    -1 3 0; #F.O
]
sinal = [
    "≤",
    "≥",
]

d = [oper, nvar, nrest, indices, sinal]

# LETRA E #
oper = "Max"
nvar = 2
nrest = 2
indices = [
    3 5 15;
    5 2 10;
    3 5 0; #F.O
]
sinal = [
    "≤",
    "≤",
]

e = [oper, nvar, nrest, indices, sinal]

# LETRA F #
oper = "Max"
nvar = 2
nrest = 2
indices = [
    1 1 1;
    2 2 4;
    3 -2 0; #F.O
]
sinal = [
    "≤",
    "≥",
]

f = [oper, nvar, nrest, indices, sinal]

# LETRA G #
oper = "Max"
nvar = 2
nrest = 2
indices = [
    2 1 10;
    1 3 9;
    1 4 0; #F.O
]
sinal = [
    "≤",
    "≤",
]

g = [oper, nvar, nrest, indices, sinal]

# LETRA H #
oper = "Max"
nvar = 2
nrest = 2
indices = [
    2 2 6;
    1 3 8;
    3 2 0; #F.O
]
sinal = [
    "≤",
    "≤",
]

h = [oper, nvar, nrest, indices, sinal]

# LETRA I #
oper = "Max"
nvar = 2
nrest = 4
indices = [
    -3 1 2;
    0 1 3;
    1 2 9;
    3 1 18
    1 1 0; #F.O
]
sinal = [
    "≤",
    "≤",
    "≤",
    "≤",
]

i = [oper, nvar, nrest, indices, sinal]

# LETRA J #
oper = "Max"
nvar = 2
nrest = 4
indices = [
    0 1 4;
    1 1 6;
    1 0 3;
    5 1 18;
    1 3 0; #F.O
]
sinal = [
    "≤",
    "≤",
    "≤",
    "≤",
]

j = [oper, nvar, nrest, indices, sinal]

# LETRA K #
oper = "Max"
nvar = 2
nrest = 5
indices = [
    1 0 3;
    0 1 2;
    1 1 6;
    1 -1 3;
    3 5 30;
    3 2 0; #F.O
]
sinal = [
    "≥",
    "≥",
    "≥",
    "≤",
    "≤",
]

k = [oper, nvar, nrest, indices, sinal]

# LETRA L #
oper = "Max"
nvar = 2
nrest = 4
indices = [
    0 1 10;
    2 5 60;
    1 1 18;
    3 1 44;
    2 1 0; #F.O
]
sinal = [
    "≤",
    "≥",
    "≤",
    "≤",
]

l = [oper, nvar, nrest, indices, sinal]

# LETRA M #
oper = "Min"
nvar = 2
nrest = 2
indices = [
    2 3 9;
    4 3 12;
    1 1 0; #F.O
]
sinal = [
    "≥",
    "≥",
]

m = [oper, nvar, nrest, indices, sinal]


# LETRA N #
oper = "Min"
nvar = 2
nrest = 3
indices = [
    5 1 10;
    2 2 12;
    1 4 12;
    3 2 0; #F.O
]
sinal = [
    "≥",
    "≥",
    "≥",
]

n = [oper, nvar, nrest, indices, sinal]

# LETRA O #
oper = "Min"
nvar = 2
nrest = 4
indices = [
    1 0 3;
    0 1 2;
    1 1 6;
    1 -1 3;    
    3 2 0; #F.O
]
sinal = [
    "≥",
    "≥",
    "≥",
    "≤",
]

o = [oper, nvar, nrest, indices, sinal]

# SOLVER #
alphabet_string = collect('A':'O')
alphabet = [a,b,c,d,e,f,g,h,i,j,k,l,m,n,o]
for x in 1:length(alphabet_string)
    println("\n============ Letra ", alphabet_string[x], " ============")
    runOptimization(alphabet[x][1], alphabet[x][2],alphabet[x][3], alphabet[x][4], alphabet[x][5])
end
