import sympy as sym 
from sympy import sin, cos
from scipy import optimize

sym.init_printing(use_unicode=False, wrap_line=True, num_columns= 1000)
def round_expr(expr, num_digits):
    return expr.xreplace({n : round(n, num_digits) for n in expr.atoms(sym.Number)})

w = sym.Symbol('w') 
x = sym.Symbol('x')
y = sym.Symbol('y')
z = sym.Symbol('z')

wt,xt,yt,zt = sym.symbols('wt xt yt zt')

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
LEN = sym.Symbol("LEN")
#quaternion of the current bone
QT = sym.Matrix(m)

# f(bone_rest_pose) = bone_pose by applying the quaternion
transform = (QT.T * arm_mat).T * T
#transform =  arm_mat.T * (QT.T * T)



order = sym.Matrix([w,x,y,z])
jacobian = transform.jacobian(order)

#sym.pretty_print(jacobian.subs([(w**2 + x**2 + y**2 + z**2, LEN)]))

sym.pretty_print(round_expr(jacobian.subs([
    (w, 0.999),
    (x, -0.044),
    (y, 0.018),
    (z, 0.014),

    (xx,1),
    (xy,0),
    (xz,0),

    (yx,0),
    (yy,0),
    (yz,1),
    
    (zx,0),
    (zy,-1),
    (zz,0),

    (hx, 0),
    (hy, 0),
    (hz, 0),

    (tx, 0),
    (ty, 1),
    (tz, 0), 
]),5));


predicted = sym.Matrix([0,0,1,1]) + jacobian * sym.Matrix([wt-w, xt-x, yt-y, zt-z])


energy = ds - jacobian * sym.Matrix([wt-w, xt-x, yt-y, zt-z])
#energy = energy.dot(energy)


ja00, ja01, ja02,ja03, ja10, ja11, ja12, ja13, ja20, ja21, ja22, ja23, ja30, ja31,ja32,ja33 = sym.symbols(
    'ja00 ja01 ja02 ja03 ja10 ja11 ja12 ja13 ja20 ja21 ja22 ja23 ja30 ja31 ja32 ja33'
)

jaco_proxy = sym.Matrix([
    [ja00, ja01, ja02, ja03],
    [ja10, ja11, ja12, ja13],
    [ja20, ja21, ja22, ja23],
    [0, 0, 0, 0]
])

LEN = sym.Symbol("LEN")

energy_proxy = ds - jaco_proxy * sym.Matrix([wt-w, xt-x, yt-y, zt-z])
#energy_proxy = energy_proxy.dot(energy_proxy)

sym.pretty_print(energy.subs([(w**2 + x**2 + y**2 + z**2, LEN)]))
#sym.pretty_print(energy_proxy)












'''

#sym.pretty_print(jacobian * sym.Matrix([wt-w, xt-x, yt-y, zt-z]))
sym.pretty_print(transform.subs([
    (w, 0.999),
    (x, -0.044),
    (y, 0.018),
    (z, 0.014),

    # (w, 1),
    # (x, 0.2),
    # (y, 0),
    # (z, 0.2),

    (xx,1),
    (xy,0),
    (xz,0),

    (yx,0),
    (yy,0),
    (yz,1),
    
    (zx,0),
    (zy,-1),
    (zz,0),

    (hx, 0),
    (hy, 0),
    (hz, 0),

    (tx, 0),
    (ty, 1),
    (tz, 0),    
]))


sym.pretty_print(predicted.subs([
    (wt, 0.999),
    (xt, -0.044),
    (yt, 0.018),
    (zt, 0.014),

    (w, 1),
    (x, 0),
    (y, 0),
    (z, 0),

    (xx,1),
    (xy,0),
    (xz,0),

    (yx,0),
    (yy,0),
    (yz,1),
    
    (zx,0),
    (zy,-1),
    (zz,0),

    (hx, 0),
    (hy, 0),
    (hz, 0),

    (tx, 0),
    (ty, 1),
    (tz, 0), 
]))


sym.pretty_print(round_expr(energy.subs([
    (w, 0.999),
    (x, -0.044),
    (y, 0.018),
    (z, 0.014),

    # (wt, 1),
    # (xt, 0),
    # (yt, 0),
    # (zt, 0),

    (dsx,0.0295),
    (dsy,-0.08736),
    (dsz,0.0043),

    (xx,1),
    (xy,0),
    (xz,0),

    (yx,0),
    (yy,0),
    (yz,1),
    
    (zx,0),
    (zy,-1),
    (zz,0),

    (hx, 0),
    (hy, 0),
    (hz, 0),

    (tx, 0),
    (ty, 1),
    (tz, 0), 
]),3), full_prec=False)

objective = round_expr(energy.subs([
    (w, 0.999),
    (x, -0.044),
    (y, 0.018),
    (z, 0.014),

    # (wt, 1),
    # (xt, 0),
    # (yt, 0),
    # (zt, 0),

    (dsx,0.0295),
    (dsy,-0.08736),
    (dsz,0.0043),

    (xx,1),
    (xy,0),
    (xz,0),

    (yx,0),
    (yy,0),
    (yz,1),
    
    (zx,0),
    (zy,-1),
    (zz,0),

    (hx, 0),
    (hy, 0),
    (hz, 0),

    (tx, 0),
    (ty, 1),
    (tz, 0), 
]),3)

sym.pretty_print(objective)


f = sym.lambdify([wt,xt,yt,zt], objective, 'numpy')

def obj(x):
    return f(x[0],x[1],x[2],x[3])

#print (optimize.minimize(obj, [0.999,-0.044,0.018,0.014]))

'''