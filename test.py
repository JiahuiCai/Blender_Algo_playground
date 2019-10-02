import sympy as sym
import cProfile


w = sym.Symbol('w') 
x = sym.Symbol('x')
y = sym.Symbol('y')
z = sym.Symbol('z')


i = sym.Symbol('i')
j = sym.Symbol('j')
k = sym.Symbol('k')
a = sym.Symbol('a')
b = sym.Symbol('b')
c = sym.Symbol('c')


sym.init_printing(use_unicode=False, wrap_line=True, num_columns= 1000)
# expr = ((w + x * i + y * j + z * k)/((w**2 + x**2 +y**2+z**2)**0.5)) * (a*i + b*j +c*k) * ((w - x * i - y * j - z * k)/((w**2 + x**2 +y**2+z**2)**0.5))

# imm_expr = expr.subs([(w, 0.854),(x, -0.146),(y, 0.354),(z, 0.354),(a, 0.0),(b, 0.0),(c, 1.0)])

# print(imm_expr)
# print("-------------------------------")
# print(sym.expand(imm_expr))


# print("-------------------------------")
# print((expr))

# print("-------------------------------")

# print(sym.expand(expr))


n = 1/(w**2 + x**2 + y**2 + z**2) 

#n = 1

m = []
m.append([])

m[0].append (1 - n * 2 * (y**2 + z**2))
m[0].append (n * 2 * (x*y + w*z))
m[0].append (n * 2 * (x*z + w*y))
m[0].append (0)

m.append([])
m[1].append (n * 2 * (x*y + w*z))
m[1].append (1 - n * 2 * (x**2 + z**2))
m[1].append (n * 2 * (y*z - w*x))
m[1].append (0)

m.append([])
m[2].append( n * 2 * (x*z - w*y))
m[2].append( n * 2 * (y*z + w*x))
m[2].append( 1 - n * 2 * (x**2 + y**2))
m[2].append (0)

m.append([])
m[3] = [0,0,0,1]

A = sym.Matrix(m)
sym.pretty_print(A)
vec = sym.Matrix([0,0,1,0])  # col vec

#sym.pretty_print(vec)

expr = A * vec


s = expr.subs([(w, 0.854),(x, -0.546),(y, 0.354),(z, 0.354)])
sp = expr.subs([(w, 0.854),(x, -0.55),(y, 0.354),(z, 0.354)])
dx = -0.55 - (-0.546)
ds = sp - s
sym.pretty_print(s)
sym.pretty_print(sp)
sym.pretty_print(ds)


point = sym.Matrix([a,b,c,0])



#sym.print_python(m[0][0])

#ax = [[sym.diff(m[i][j],x) for i in range(4)] for j in range(4)]

print("------------------------------")
print("------------------------------")
print("------------------------------")
print("------------------------------")
print("------------------------------")

sym.pretty_print(sym.diff(m[2][2],x))
print("------------------------------")
print("------------------------------")
print("------------------------------")
print("------------------------------")
print("------------------------------")
ax = sym.simplify(sym.diff(A,x))
sym.pretty_print(ax)
print("------------------------------")
print("------------------------------")
print("------------------------------")
print("------------------------------")
print("------------------------------")

newA = A + dx * ax 

expr_appr = newA * vec


sp_appr = expr_appr.subs([(w, 0.854),(x, -0.546),(y, 0.354),(z, 0.354)])
#sym.pretty_print(sp_appr)
#sym.pretty_print(sp)


expr_appr = A * vec + (-0.005) * ax * vec
sp_appr = expr_appr.subs([(w, 0.854),(x, -0.546),(y, 0.354),(z, 0.354)]) 

#sym.pretty_print(sp_appr)
#sym.pretty_print(sp)


xp = sym.Symbol('xp')
energy = ((ds - (xp - (x)) * ax * vec))
energy = energy.dot(energy)


sym.pretty_print(energy)
print("\nwtf\n")
sym.pretty_print(sym.collect(sym.expand(energy), xp).coeff(xp,3).simplify())
print("\nwtf3\n")
sym.pretty_print(sym.collect(sym.expand(energy), xp).coeff(xp,2).simplify())
print("\nwtf2\n")
sym.pretty_print(sym.collect(sym.expand(energy), xp).coeff(xp,1).simplify())
print("\nwtf1\n")
sym.pretty_print(sym.collect(sym.expand(energy), xp).coeff(xp,0).simplify())
print("\nwtf0\n")
#cProfile.run('sym.solveset(sym.Eq(energy.subs([(w, 0.854),(x, -0.546),(y, 0.354),(z, 0.354)]),4.17767006425792e-10), xp, domain = sym.S.Reals)')

sym.pretty_print(sym.simplify(energy.subs([(w, 0.854),(x, -0.546),(y, 0.354),(z, 0.354)])))

sym.pretty_print(sym.solveset(sym.Eq(energy.subs([(w, 0.854),(x, -0.546),(y, 0.354),(z, 0.354)]),4.17767006425792e-10), xp, domain = sym.S.Reals))



zero_ener = energy.subs(([(w, 0.854),(x, -0.546),(y, 0.354),(z, 0.354), (xp, -0.55)]))

sym.pretty_print(zero_ener)


