import sympy as sym 

sym.init_printing(use_unicode=False, wrap_line=True, num_columns= 1000)
w = sym.Symbol('w') 
x = sym.Symbol('x')
y = sym.Symbol('y')
z = sym.Symbol('z')

xt = sym.Symbol('xt')

hx = sym.Symbol('hx')
hy = sym.Symbol('hy')
hz = sym.Symbol('hz')
# these are the arm_space_mat
xx,xy,xz,yx,yy,yz,zx,zy,zz = sym.symbols('xx xy xz yx yy yz zx zy zz')

arm_mat = sym.Matrix([
    [xx,xy,xz, 0],
    [yx,yy,yz, 0],
    [zx,zy,zz, 0],
    [hx,hy,hz,1]
])


# these are bone space shit
tx = sym.Symbol('tx')
ty = sym.Symbol('ty')
tz = sym.Symbol('tz')
T = sym.Matrix([tx,ty,tz,1])

Ori = sym.Matrix([0,0,0,1])


dsx = sym.Symbol('dsx')
dsy = sym.Symbol('dsy')
dsz = sym.Symbol('dsz')
ds = sym.Matrix([dsx, dsy, dsz, 0])

# working with full 4D input space, not just unit quaternions
n = 1/(w**2 + x**2 + y**2 + z**2) 
#n = 1

m = []
m.append([])

m[0].append (1 - n * 2 * (y**2 + z**2))
m[0].append (n * 2 * (x*y - w*z))
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

#quaternion of the current bone
QT = sym.Matrix(m)

# f(bone_rest_pose) = bone_pose by applying the quaternion
transform = (QT.T * arm_mat).T * T

# partial derivatives of the Quaternion matrix with respect to the x component
dQTdx = sym.simplify(sym.diff(QT.T,x))
sym.pretty_print(dQTdx)

# using partial derivatives to estimate the output of a new x compoenent xt T(xt) = T(x) + dx * (dQT/dx) * arm_space_bone_mat * bone_space_tail
predicted = transform + (xt - x) * ((dQTdx * arm_mat).T * T)

# ds is the delta vec between current location and the target location of the bone tail
# minimize this energy gives a best xt that produces a quaternion that rotates the bone to the target location
# a quadratic energy, solvable with many solvers
J = ((dQTdx * arm_mat).T * T ) 

sym.pretty_print(J)
print("ccode")
sym.printing.print_ccode(J[0])
print("ccode---end")
Jx, Jy, Jz = sym.symbols('Jx, Jy, Jz')
J_proxy = sym.Matrix([Jx,Jy,Jz, 0])

energy_proxy = ds - (xt - x) * J_proxy
energy_proxy_x = energy_proxy[0] ** 2
energy_proxy = energy_proxy.dot(energy_proxy)
sym.pretty_print(sym.collect(energy_proxy.expand(),xt))


print("constant")
sym.printing.print_ccode(sym.collect(energy_proxy.expand(),xt).coeff(xt, 0))
print("linear")
sym.printing.print_ccode(sym.collect(energy_proxy.expand(),xt).coeff(xt, 1))
print("quadratic")
sym.printing.print_ccode(sym.collect(energy_proxy.expand(),xt).coeff(xt, 2))


energy = ds - (xt - x) * J
energy = energy.dot(energy)

LEN = sym.Symbol('LEN')

#sym.pretty_print (sym.collect(energy.subs([(w**2 + x**2 + y**2 + z**2, LEN)]).expand(), xt).coeff(xt, 0).simplify()) 
#sym.pretty_print ((energy.subs([(w**2 + x**2 + y**2 + z**2, LEN)]))) 
#sym.printing.print_ccode((energy.subs([(w**2 + x**2 + y**2 + z**2, LEN)])))
########### test cases ######################

'''
sym.pretty_print ( transform.subs ([
    (w, 0.707),
    (x, 0.707),
    (y, 0.000),
    (z, 0.000),
    (xx, 0),
    (xy, 1),
    (xz, 0),
    (yx, -0.707),
    (yy, 0),
    (yz, 0.707),
    (zx, 0.707),
    (zy, 0 ),
    (zz, 0.707),
    (hx, -1),
    (hy, 0),
    (hz, 1),
    (tx, 0),
    (ty, 1.414),
    (tz, 0),
])) 


sym.pretty_print(predicted.subs([
    (w, 0.707),
    (x, 0.707),
    (y, 0.000),
    (z, 0.000),
    (xx, 0),
    (xy, 1),
    (xz, 0),
    (yx, -0.707),
    (yy, 0),
    (yz, 0.707),
    (zx, 0.707),
    (zy, 0 ),
    (zz, 0.707),
    (hx, -1),
    (hy, 0),
    (hz, 1),
    (tx, 0),
    (ty, 1.414),
    (tz, 0),
    (xt, 0.75)
]))

sym.pretty_print(energy.subs([
    (w, 0.707),
    (x, 0.707),
    (y, 0.000),
    (z, 0.000),
    (xx, 0),
    (xy, 1),
    (xz, 0),
    (yx, -0.707),
    (yy, 0),
    (yz, 0.707),
    (zx, 0.707),
    (zy, 0 ),
    (zz, 0.707),
    (hx, -1),
    (hy, 0),
    (hz, 1),
    (tx, 0),
    (ty, 1.414),
    (tz, 0),
    (xt, 0.75),
    (dsx, 0.057),
    (dsy, 0),
    (dsz, -0.061)
]))


sym.pretty_print(sym.solveset(sym.Eq(energy.subs([
    (w, 0.707),
    (x, 0.707),
    (y, 0.000),
    (z, 0.000),
    (xx, 0),
    (xy, 1),
    (xz, 0),
    (yx, -0.707),
    (yy, 0),
    (yz, 0.707),
    (zx, 0.707),
    (zy, 0 ),
    (zz, 0.707),
    (hx, -1),
    (hy, 0),
    (hz, 1),
    (tx, 0),
    (ty, 1.414),
    (tz, 0),
    #(xt, 0.75),
    (dsx, 0.057),
    (dsy, 0),
    (dsz, -0.061)
]), 0.0001), xt, domain = sym.S.Reals))



sym.pretty_print ( J.subs ([
    (w, 0.707),
    (x, 0.707),
    (y, 0.000),
    (z, 0.000),
    (xx, 0),
    (xy, 1),
    (xz, 0),
    (yx, -0.707),
    (yy, 0),
    (yz, 0.707),
    (zx, 0.707),
    (zy, 0 ),
    (zz, 0.707),
    (hx, -1),
    (hy, 0),
    (hz, 1),
    (tx, 0),
    (ty, 1.414),
    (tz, 0),
])) 


'''

sym.pretty_print ( energy_proxy.subs ([
    (x, 0.707),
    (xt, 0.75),
    (dsx, 0.057),
    (dsy, 0),
    (dsz, -0.061),
    (Jx,1.414),
    (Jy,0),
    (Jz,-1.414), 
])) 

sym.pretty_print(sym.solveset(sym.Eq(energy_proxy.subs([
    (x, 0.707),
    #(xt, 0.75),
    (dsx, 0.057),
    (dsy, 0),
    (dsz, -0.061),
    (Jx,1.414),
    (Jy,0),
    (Jz,-1.414), 
]), 0.0001), xt, domain = sym.S.Reals))


sym.pretty_print(sym.solveset(sym.Eq(energy_proxy_x.subs([
    (x, 0.707),
    #(xt, 0.75),
    (dsx, 0.057),
    (Jx,1.414),
]), 0.0001), xt, domain = sym.S.Reals))