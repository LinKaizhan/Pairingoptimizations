ells = [ 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
	103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211,
	223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347,
	349, 353, 359, 367, 373, 587];
p = 4 * prod(ells) - 1;
G.<a>=GF(p^2,name='a',modulus=x^2+1);
F = GF(p);
A = 6;
#A = 3660531689687144537798649124765796265922835935302812894295885379832752438089867533337299606759071662303728117187051820332336016283001207986531180645517342;
#A = 0;
E = EllipticCurve(G,[0,A,0,1,0]);
N = 0x10A59877399FC365649BD0A396FF5CCF84EAAEB92E2E947EA2113B6E1D850F5F7C67;
COFACTOR = ceil((p+1)/N);
ellsN = ells[41:74];
u = G(2);
sqrt_u = sqrt(u);


# Point doubling and point addition
def xDBLADD(X1, Z1, X2, Z2, X3, Z3, A24):
    t0 = X1 + Z1; t1 = X1 - Z1;     #t0 = XP+ZP, t1 = XP-ZP
    X4 = t0**_sage_const_2 ;        #X4 = (XP+ZP)^2
    t2 = X2 - Z2;                   #t2 = XQ-ZQ
    X5 = X2 + Z2;                   #X5 = XQ+ZQ
    t0 = t0 * t2;                   #t0 = (XP+ZP)*(XQ-ZQ)
    Z4 = t1**_sage_const_2 ;        #Z4 = (XP-ZP)^2
    t1 = t1 * X5; t2 = X4 -Z4;      #t1 = (XP-ZP)*(XQ+ZQ), t2 = 4*XP*ZP
    X4 = X4 * Z4;                   #X4 = (XP+ZP)^2*(XP-ZP)^2
    X5 = A24 * t2;                  #X5 = A24 * 4*XP*ZP
    Z5 = t0 - t1;                   #Z5 = 2*(XQZP-XPZQ)
    Z4 = X5 + Z4;                   #Z4 = A24 * 4*XP*ZP+(XP-ZP)^2
    X5 = t0 + t1;                   #X5 = 2*(XPZP-ZPZQ)
    Z4 = Z4 * t2;                   #Z4 = 4*XP*ZP*(A24 * 4*XP*ZP+(XP-ZP)^2)
    Z5 = Z5**_sage_const_2;         #Z5 = 4*(XQZP-XPZQ)^2
    X5 = X5**_sage_const_2;         #X5 = 4*(XPZP-ZPZQ)^2
    Z5 = X3 * Z5;
    X5 = Z3 * X5;
    if Z3 == 1:
        cost_count[2] = cost_count[2] + 6;
        cost_count[3] = cost_count[3] + 4;
    else:
        cost_count[2] = cost_count[2] + 7;
        cost_count[3] = cost_count[3] + 4;
    return X4, Z4, X5, Z5;

# Montgomery ladder
def Montgomery_ladder(ell, A, XP, ZP):
    X2 = XP;
    Z2 = ZP;
    A24 = A + G(2);
    C24 = G(4);
    X1 = 1; Z1 = 0;
    ellbits = bin(ell);
    ellbits = tuple(ellbits);
    ellbits = ellbits[2:];
    l = len(ellbits);
    A24 = A24/C24;
    for j in range(l):
        if ellbits[j] == '1':
            [X2, Z2, X1, Z1] = xDBLADD(X2, Z2, X1, Z1, XP, ZP, A24);
        else:
            [X1, Z1, X2, Z2] = xDBLADD(X1, Z1, X2, Z2, XP, ZP, A24);
    return X1, Z1, X2, Z2;


# Batch cofactor multiplication
# Output:{Ui},the order of each point is ln,ln-1,...,l1.
def batch_cof_mul(PXi,PZi,XP,ZP,ind,n,A,L):
    if n == 1:
        PXi.append(XP);
        PZi.append(ZP);
        return 1;
    m=floor(n/2);
    u=1;
    v=1;
    for i in range(m):
        u=L[ind[i]]*u;
    for j in range(m,n):
        v=L[ind[j]]*v;
    R0X, R0Z, R1X, R1Z = Montgomery_ladder(u, A, XP, ZP);
    batch_cof_mul(PXi,PZi,R0X,R0Z,[ind[i] for i in range(m,n)],n-m,A,L);
    R0X, R0Z, R1X, R1Z = Montgomery_ladder(v, A, XP, ZP);
    batch_cof_mul(PXi,PZi,R0X,R0Z,[ind[i] for i in range(m)],m,A,L);

def recover_y3(X0, Y0, Z0, X1, Z1, X2, Z2, A):
    x = X0;
    y = Y0;
    z = Z0;
    t1 = x * Z1;
    t2 = X1 * z;
    t3 = t2 - t1;           #t3 = X1*z-x*Z1
    t2 = t2 + t1;           #t2 = X1*z+x*Z1
    t3 = t3^2;              #t3 = (X1*z-x*Z1)^2
    t3 = t3 * X2;           #t3 = (X1*z-x*Z1)^2*X2
    t1 = 2 * Z1;            #t1 = 2*Z1
    t1 = t1 * A;            #t1 = 2*Z1*A
    t1 = t1 * z;            #t1 = 2*Z1*A*z
    t2 = t1 + t2;           #t2 = 2*Z1*A*z+X1*z+x*Z1
    t4 = x * X1;            #t4 = x*X1
    t4 = t4 + Z1*z;         #t4 = x*X1+Z1*z
    t2 = t2 * t4;           #t2 = (2*Z1*A*z+X1*z+x*Z1)*(x*X1+Z1*z)
    t1 = t1 * z;            #t1 = 2*Z1*A*z^2
    t1 = t1 * Z1;           #t1 = 2*Z1*A*z^2*Z1
    t2 = t2 - t1;           #t2 = (2*Z1*A*z+X1*z+x*Z1)*(x*X1+Z1*z)-2*Z1*A*z^2*Z1
    t2 = t2 * Z2;           #t2 = ((2*Z1*A*z+X1*z+x*Z1)*(x*X1+Z1*z)-2*Z1*A*z^2*Z1)*Z2
    Y = t2 - t3;
    t1 = 2 * y * z;
    t1 = t1 * Z1;
    t1 = t1 * Z2;
    X = X1 * t1;
    Z = Z1 * t1;
    if X == 0 and Y ==0 and Z == 0:
        Y = G(1);
    if Z0 == 1:
        cost_count[2] = cost_count[2] + 11;
    else:
        cost_count[2] = cost_count[2] + 16;
    cost_count[3] = cost_count[3] + 1;
    return X, Y, Z;

# elligator: randomly generate a point
def elligator(A, E):
  flag_2 = 0;
  while flag_2 == 0:
    r = randint(0,p);
    r = F(r);
    tmp = F(u) * r^2;
    cost_count[2] = cost_count[2] + 1;
    cost_count[3] = cost_count[3] + 1;
    cond1 = 1 + tmp;
    cond2l = A^2 * tmp;
    cost_count[2] = cost_count[2] + 1;
    cost_count[3] = cost_count[3] + 1;
    cond2r = cond1^2;
    cost_count[3] = cost_count[3] + 1;
    if cond1 != 0 and cond2l != cond2r:
        v = -A/cond1;
        cost_count[4] = cost_count[4] + 1;
        v = F(v);
        e = is_square(v * (v^2 + A*v + F(1)));
        # In fact we do not need to compute e as y1 should be computed
        x1 = G(v); x2 = -v - A;
        y1 = sqrt(x1 * (x1^2 + A * x1 + 1));
        cost_count[2] = cost_count[2] + 2;
        cost_count[3] = cost_count[3] + 1;
        cost_count[2] = cost_count[2] + 255;
        cost_count[3] = cost_count[3] + 509;
        y2 = y1*sqrt_u*r;
        # We can select a small r to save one multiplication
        cost_count[2] = cost_count[2] + 1;
        if e == True:
            return E([x1,y1]), E([x2,y2]);
        else:
            return E([x2,y2]), E([x1,y1]);
        flag_2 = 1;

# Translate a string to its naf form
def Trans(N):
    N_str = N.str(base=2);
    l = len(N_str);
    string_naf = [];
    for j in range(l+1):
        string_naf.append(0);
    j = 0;
    while(j < l):
        if N_str[l-j-1] == '1':
            k = j;
            while(k != l-1 and N_str[l-k-2] == '1'):
                k = k + 1;
            if k != j:
                string_naf[l-j] = - 1;
                string_naf[l-k-1] = 1;
                j = k;
            else:
                string_naf[l-j] = 1;
        j = j + 1;
    j = 0;
    while(j < l):
        if string_naf[l-j] == 1:
            k = j;
            while(k != l and string_naf[l-k-1] == 1):
                k = k + 1;
            if k != j:
                string_naf[l-j] = - 1;
                for jj in range(j,k):
                    string_naf[l-jj-1] = 0;
                string_naf[l-k-1] = 1;
                j = k;
            if string_naf[l-j-1] == -1:
                string_naf[l-j-1] = 0;
                string_naf[l-j] = -1;
                j = j - 1;
        j = j + 1;
    if string_naf[0] == 0:
        for j in range(l):
            string_naf[j] = string_naf[j+1];
        string_naf.pop();
    return string_naf;

def Lucassequences(ft,ell):
	v0 = 2;
	tmp1 = ft;
	v1 = tmp1;
	tmp2 = G(2);
	string = ell.str(base=2);
	Len = len(string);
	for j in range(Len):
		if string[j] == '1':
			v0 = v0 * v1 - tmp1; cost_count[2] = cost_count[2] + 1;
			v1 = v1 * v1 - tmp2; cost_count[3] = cost_count[3] + 1;
		else:
			v1 = v0 * v1 - tmp1; cost_count[3] = cost_count[3] + 1;
			v0 = v0 * v0 - tmp2; cost_count[2] = cost_count[2] + 1;
	return v0;

def ADD_Projective(X1,Y1,Z1,X2,Y2,Z2,AA,BB):
    b3 = BB*3;
    t0 = X1*X2;
    t1 = Y1*Y2;
    t2 = Z1*Z2;
    t3 = X1+Y1;
    t4 = X2+Y2;
    t3 = t3*t4;
    t4 = t0+t1;
    t3 = t3-t4;
    t4 = X1+Z1;
    t5 = X2+Z2;
    t4 = t4*t5;
    t5 = t0+t2;
    t4 = t4-t5;
    t5 = Y1+Z1;
    X3 = Y2+Z2;
    t5 = t5*X3;
    X3 = t1+t2;
    t5 = t5-X3;
    Z3 = AA*t4;
    X3 = b3*t2;
    Z3 = X3+Z3;
    X3 = t1-Z3;
    Z3 = t1+Z3;
    Y3 = X3*Z3;
    t1 = t0+t0;
    t1 = t1+t0;
    t2 = AA*t2;
    t4 = b3*t4;
    t1 = t1+t2;
    t2 = t0-t2;
    t2 = AA*t2;
    t4 = t4+t2;
    t0 = t1*t4;
    Y3 = Y3+t0;
    t0 = t5*t4;
    X3 = t3*X3;
    X3 = X3-t0;
    t0 = t3*t1;
    Z3 = t5*Z3;
    Z3 = Z3+t0;
    if Z2 == 1 or Z1 == 1:
        cost_count[2] = cost_count[2] + 16;
    else:
        cost_count[2] = cost_count[2] + 17;
    return X3,Y3,Z3;

# point doubling in modified Jacobian coordinates
# Cost:3m+5s
def DBL(xx,yy,zz,tt,AA):
    t1 = xx^2;
    t2 = 2*yy^2;
    t3 = t2^2;
    t4 = 2*t3;
    t5 = (xx+t2)^2-t1-t3;
    lam = 3*t1+tt;
    x3 = lam^2-2*t5;
    y3 = lam*(t5-x3)-t4;
    z3 = 2*yy*zz;
    t3 = 2*t4*tt;
    cost_count[2] = cost_count[2] + 3;
    cost_count[3] = cost_count[3] + 5;
    return x3,y3,z3,t3,lam;

# point doubling in Jacobian coordinates
# Cost: 2m+8s (2m+7s when Z12=Z1^2 is given)
def DBL_Jacobian(X1,Y1,Z1,Z12,AA):
    XX = X1^2;      # 1s
    YY = Y1^2;      # 1s
    YYYY = YY^2;    # 1s
    ZZ = Z12;
    t0 = X1+YY;
    t1 = t0^2;      # 1s
    t2 = t1-XX;
    t3 = t2-YYYY;
    S = 2*t3;
    t4 = ZZ^2;      # 1s
    t5 = AA*t4;     # 1m
    t6 = 3*XX;
    M = t6+t5;
    t7 = M^2;       # 1s
    t8 = 2*S;
    T = t7-t8;
    X3 = T;
    t9 = S-T;
    t10 = 8*YYYY;
    t11 = M*t9;     # 1m
    Y3 = t11-t10;
    t12 = Y1+Z1;
    t13 = t12^2;    # 1s
    t14 = t13-YY;
    Z3 = t14-ZZ;
    cost_count[2] = cost_count[2] + 2;
    cost_count[3] = cost_count[3] + 7;
    return X3,Y3,Z3,M;

# point addition in modified Jacobian coordinates (Z2==1)
# Cost: 8m+6s
def ADD(xx,yy,zz,x1,y1,AA):
    t1 = zz^2;              # 1s
    t2 = x1*t1-xx;          # 1m
    t3 = t2^2;              # 1s
    t4 = 4*t3;
    t5 = t2*t4;             # 1m
    t6 = y1*zz*t1-yy;       # 2m
    t7 = 2 * t6;
    t8 = t4*xx;             # 1m
    x3 = t7^2-t5-2*t8;      # 1s
    y3 = t7*(t8-x3)-2*yy*t5;# 2m
    z3 = (zz+t2)^2-t1-t3;   # 1s
    t3 = AA*z3^4;           # 2s+1m
    cost_count[2] = cost_count[2] + 8;
    cost_count[3] = cost_count[3] + 6;
    return x3,y3,z3,t3,t1,t2,t6;

# point addition in Jacobian coordinates (Z2==1)
# Cost: 7m+4s
def ADD_Jacobian(X1,Y1,Z1,X2,Y2):
    Z1Z1 = Z1^2;    # 1s
    U2 = X2*Z1Z1;   # 1m
    t0 = Z1*Z1Z1;   # 1m
    S2 = Y2*t0;     # 1m
    H = U2-X1;
    HH = H^2;       # 1s
    I = 4*HH;
    J = H*I;        # 1m
    t1 = S2-Y1;
    r = 2*t1;
    V = X1*I;       # 1m
    t2 = r^2;       # 1s
    t3 = 2*V;
    t4 = t2-J;
    X3 = t4-t3;
    t5 = V-X3;
    t6 = Y1*J;      # 1m
    t7 = 2*t6;
    t8 = r*t5;      # 1m
    Y3 = t8-t7;
    t9 = Z1+H;
    t10 = t9^2;     # 1s
    t11 = t10-Z1Z1;
    Z3 = t11-HH;
    cost_count[2] = cost_count[2] + 7;
    cost_count[3] = cost_count[3] + 4;
    return X3,Y3,Z3,Z1Z1,t0,H,t1;

# Miller line evaluation for a doubling step
# Cost: 1M+1S+4m+1s
def line_fun_double(xx,yy,zz,X,Y,tmp,f):
    flag = 0;
    t1 = zz^2;                  # 1s
    t2 = t1 * zz;               # 1m
    t3 = t1 * X;                # 1m
    t3 = t3 - xx;
    t3 = t3 * tmp;              # 1m
    t2 = t2 * Y;                # 1m
    t3 = t3 - t2;
    t3 = t3 - yy;
    f = f^2 * t3;               # 1M+1S
    if flag == 1:
        f = 1;
    cost_count[0] = cost_count[0] + 1;
    cost_count[1] = cost_count[1] + 1;
    cost_count[2] = cost_count[2] + 4;
    cost_count[3] = cost_count[3] + 1;
    return f,t1;

# Miller line evaluation for an addition step
# Cost: 1M+3m
def line_fun_addition(xx,yy,zz,x1,y1,tmp2,tmp6,X,Y,f):
    t1 = X-x1;
    t1 = tmp6 * t1;     # 1m
    t2 = Y - y1;
    t2 = tmp2 * t2;     # 2m
    t1 = t1 - t2;
	#t1 = (z^4*y1+y*z)*(X-x1)-z^2*(z^2*x1-x)*(Y-y1);
    if t2 != 0:
        # ((z^4*y1+y*z)*(X-x1)-z^2*(z^2*x1-x)*(Y-y1))/((z^2*X-x)*(z^2*x1-x))
        f = f * t1;     # 1M
    else:
        f = 1;
    cost_count[0] = cost_count[0] + 1;
    cost_count[2] = cost_count[2] + 3;
    return f;

# Miller line evaluation for a double-and-addition step in modified Jacobian coordinates
# Cost: 2M+1S+6m
def line_fun_dbladd(xx,yy,zz,x1,y1,tmp,tmp1,tmp2,tmp6,X,Y,f):
    #t1 = (tmp*(X*zz^2-xx)-(Y*zz^3+yy))*tmp2;
    #t2 = X-x1;
    #t2 = (tmp6)*(X*zz^2-xx)+(Y*zz^3+yy)*tmp2;
    t1 = tmp1;
    t2 = t1 * zz;                   # 1m
    t3 = t1 * X - xx;               # 1m
    t4 = tmp * t3;                  # 1m
    t5 = Y * t2;                    # 1m
    t6 = t5 + yy;
    t7 = t4 - t6;
    t8 = t7;
    t9 = 1;
    t10 = tmp6 * t3;                # 1m
    t11 = t6 * tmp2;                # 1m
    t12 = t10 + t11;
    if t12 != 0:
        t1 = t8*(t12.conjugate());  # 1M
        f = f^2 * t1;               # 1S+1M
    else:
        f = 1;
    cost_count[0] = cost_count[0] + 2;
    cost_count[1] = cost_count[1] + 1;
    cost_count[2] = cost_count[2] + 6;
    return f;

# Miller line evaluation for a double-and-subtraction step in modified Jacobian coordinates
# Cost: 2M+1S+6m
def line_fun_dblmin(xx,yy,zz,x1,y1,tmp,tmp1,tmp2,tmp6,X,Y,f):
    #t1 = (tmp*(X*zz^2-xx)-(Y*zz^3+yy))*tmp2;
    #t2 = (tmp6)*(X*zz^2-xx)+(Y*zz^3+yy)*tmp2;
    t1 = tmp1;
    t2 = t1 * zz;               # 1m
    t3 = t1 * X - xx;           # 1m
    t4 = tmp * t3;              # 1m
    t5 = Y * t2;                # 1m
    t6 = t5 + yy;
    t7 = t4 - t6;
    t8 = t7;
    t9 = tmp6 * t3;             # 1m
    t10 = t6 * tmp2;            # 1m
    t11 = t9 + t10;
    if t11 != 0:
        t1 = t8*t11.conjugate();# 1M
        f = f^2 * t1;           # 1S+1M
    else:
        f = 1;
    cost_count[0] = cost_count[0] + 2;
    cost_count[1] = cost_count[1] + 1;
    cost_count[2] = cost_count[2] + 6;
    return f;

# Miller line evaluation for a double-and-addition step in Jacobian coordinates
# Cost: 2M+1S+6m
def line_fun_Jacobian_dbladd(xx,yy,zz,x1,y1,tmp,Z1Z1,t0,H,r,X,Y,f):
    #t1 = (tmp*(X*zz^2-xx)-(Y*zz^3+yy))*tmp2;
    #t2 = X-x1;
    #t3 = (tmp6)*(X*zz^2-xx)+(Y*zz^3+yy)*tmp2;
    t1 = Z1Z1;
    t2 = t0;
    t3 = t1 * X - xx;               # 1m
    t4 = tmp * t3;                  # 1m
    t5 = Y * t2;                    # 1m
    t6 = t5 + yy;
    t7 = t4 - t6;
    t8 = t7;
    t9 = r * t3;                    # 1m
    t10 = t6 * H;                   # 2m
    t11 = t9 + t10;
    if t11 != 0:
        t1 = t8*(t11.conjugate());  # 1M
        f = f^2 * t1;               # 1M+1S
    else:
        f = 1;
    cost_count[0] = cost_count[0] + 2;
    cost_count[1] = cost_count[1] + 1;
    cost_count[2] = cost_count[2] + 6;
    return f;

# Miller line evaluation for a double-and-subtraction step in Jacobian coordinates
#Cost: 2M+1S+6m
def line_fun_Jacobian_dblmin(xx,yy,zz,x1,y1,tmp,Z1Z1,t0,H,r,X,Y,f):
    #t1 = (tmp*(X*zz^2-xx)-(Y*zz^3+yy))*tmp2;
    #t2 = (tmp6)*(X*zz^2-xx)+(Y*zz^3+yy)*tmp2;
    t1 = Z1Z1;
    t2 = t0;
    t3 = t1 * X - xx;           # 1m
    t4 = tmp * t3;              # 1m
    t5 = Y * t2;                # 1m
    t6 = t5 + yy;
    t7 = t4 - t6;
    t8 = t7;
    t9 = r * t3;                # 1m
    t10 = t6 * H;               # 2m
    t11 = t9 + t10;
    if t11 != 0:
        t1 = t8*t11.conjugate();# 1M
        f = f^2 * t1;           # 1M+1S
    else:
        f = 1;
    cost_count[0] = cost_count[0] + 2;
    cost_count[1] = cost_count[1] + 1;
    cost_count[2] = cost_count[2] + 6;
    return f;

# Miller line evaluation for a quadrupling step
# Cost: 2M+2S+7m+1s
def line_fun_qua(xx2,yy2,zz2,X,Y,tmp1,tmp2,f):
    flag = 0;
    t1 = zz2^2;                 # 1s
    t2 = t1 * zz2;              # 1m
    t3 = tmp2;
    t4 = yy2;
    f = f^2;                    # 1S
    t5 = t1 * X;                # 1m
    t5 = t5 - xx2;
    t6 = t5 * tmp1;             # 1m
    t7 = Y * t2;                # 1m
    t7 = t7 + yy2;
    t6 = t6 - t7;
    f = f * t6;                 # 1M
    f = f^2;                    # 1S
    t5 = t5 * t3;               # 1m
    t6 = t4 * t7;               # 2m
    t6 = t6 + t6;
    t5 = t5 + t6;
    if t5 != 0:
        f = f*t5.conjugate();   # 1M
    else:
        flag = 1;
    if flag == 1:
        f = 1;
    cost_count[0] = cost_count[0] + 2;
    cost_count[1] = cost_count[1] + 2;
    cost_count[2] = cost_count[2] + 7;
    cost_count[3] = cost_count[3] + 1;
    return f;

# Miller line evaluation for a quadrupling step
# Cost: 2M+2S+7m
def line_fun_Jacobian_qua(xx2,yy2,zz2,zzs,X,Y,tmp1,tmp2,f):
    flag = 0;
    t1 = zzs;
    t2 = t1 * zz2;              # 1m
    t3 = tmp2;
    t4 = yy2;
    f = f^2;                    # 1S
    t5 = t1 * X;                # 1m
    t5 = t5 - xx2;
    t6 = t5 * tmp1;             # 1m
    t7 = Y * t2;                # 1m
    t7 = t7 + yy2;
    t6 = t6 - t7;
    f = f * t6;                 # 1M
    f = f^2;                    # 1S
    t5 = t5 * t3;               # 1m
    t6 = t4 * t7;               # 2m
    t6 = t6 + t6;
    t5 = t5 + t6;
    f = f*t5.conjugate();       # 1M
    cost_count[0] = cost_count[0] + 2;
    cost_count[1] = cost_count[1] + 2;
    cost_count[2] = cost_count[2] + 7;
    return f;

def miller_modified_odd_naf_trace(X,Y,xx,yy,zz,tt,AA,exp,fin):
    f=G(1);
    string = Trans(exp);
    l=len(string);
    flag=0;
    # Note that zz == 1
    x1 = xx;
    y1 = yy;
    i = 1;
    while (i != l):
        [xx,yy,zz,tt,tmp] = DBL(xx,yy,zz,tt,AA);
        if string[i] == 1 and i != l-1:
            if xx != x1*zz^2:
                [xx2,yy2,zz2,tt,tmp1,tmp2,tmp6] = ADD(xx,yy,zz,x1,y1,AA);
                f = line_fun_dbladd(xx,yy,zz,x1,y1,tmp,tmp1,tmp2,tmp6,X,Y,f);
                xx = xx2;
                yy = yy2;
                zz = zz2;
                i = i + 1;
            else:
                flag = 1;
        else:
            if string[i] == -1 and i != l-1:
                if xx != x1*zz^2:
                    [xx2,yy2,zz2,tt,tmp1,tmp2,tmp6] = ADD(xx,yy,zz,x1,-y1,AA);
                    f = line_fun_dblmin(xx,yy,zz,x1,y1,tmp,tmp1,tmp2,tmp6,X,Y,f);
                    xx = xx2;
                    yy = yy2;
                    zz = zz2;
                    i = i + 1;
            else:
                if string[i] == 0 and string[i+1] == 0 and i <= l-2:
                    [xx2,yy2,zz2,tt,tmp2] = DBL(xx,yy,zz,tt,AA);
                    f = line_fun_qua(xx,yy,zz,X,Y,tmp,tmp2,f);
                    i = i + 2;
                    xx = xx2;
                    yy = yy2;
                    zz = zz2;
                else:
                    [f,zzs] = line_fun_double(xx,yy,zz,X,Y,tmp,f);
                    i = i + 1;
    f = f^(p-1);
    cost_count[2] = cost_count[2] + 3;
    cost_count[3] = cost_count[3] + 2;
    cost_count[4] = cost_count[4] + 1;
    f = f + f.conjugate();
    f = Lucassequences(f,ceil(p+1)/fin);
    if flag == 1:
        f = 2;
    return f;

# Pairing (using new formulas and Jacobian coordinates)
# the order of the pairing is odd (the string has naf form)
def miller_Jacobian_odd_naf_trace(X,Y,xx,yy,zz,AA,exp,fin):
    f=G(1);
    string = Trans(exp);
    l=len(string);
    flag=0;
    # Note that zz == 1
    zzs = 1;
    x1 = xx;
    y1 = yy;
    i = 1;
    while (i != l):
        [xx,yy,zz,tmp] = DBL_Jacobian(xx,yy,zz,zzs,AA);
        if string[i] == 1 and i != l-1:
            if xx != x1*zz^2:
                [xx2,yy2,zz2,Z1Z1,t0,H,r] = ADD_Jacobian(xx,yy,zz,x1,y1);
                f = line_fun_Jacobian_dbladd(xx,yy,zz,x1,y1,tmp,Z1Z1,t0,H,r,X,Y,f);
                xx = xx2;
                yy = yy2;
                zz = zz2;
                zzs = zz^2;
                cost_count[3] = cost_count[3] + 1;
                i = i + 1;
            else:
                flag = 1;
        else:
            if string[i] == -1 and i != l-1:
                if xx != x1*zz^2:
                    [xx2,yy2,zz2,Z1Z1,t0,H,r] = ADD_Jacobian(xx,yy,zz,x1,-y1);
                    f = line_fun_Jacobian_dblmin(xx,yy,zz,x1,y1,tmp,Z1Z1,t0,H,r,X,Y,f);
                    xx = xx2;
                    yy = yy2;
                    zz = zz2;
                    zzs = zz^2;
                    cost_count[3] = cost_count[3] + 1;
                    i = i + 1;
            else:
                if string[i] == 0 and string[i+1] == 0 and i < l-2:
                    zzs = zz^2;
                    cost_count[3] = cost_count[3] + 1;
                    [xx2,yy2,zz2,tmp2] = DBL_Jacobian(xx,yy,zz,zzs,AA);
                    f = line_fun_Jacobian_qua(xx,yy,zz,zzs,X,Y,tmp,tmp2,f);
                    i = i + 2;
                    xx = xx2;
                    yy = yy2;
                    zz = zz2;
                    zzs = zz2^2;
                    cost_count[3] = cost_count[3] + 1;
                else:
                    [f,zzs] = line_fun_double(xx,yy,zz,X,Y,tmp,f);
                    i = i + 1;
    if string[l-1] == 1:
        #if X != x1:
            #f = f * (X-x1);
            #cost_count[2] = cost_count[2] + 2;
        if X == x1:
            flag = 1;
    f = f^(p-1);
    cost_count[2] = cost_count[2] + 3;
    cost_count[3] = cost_count[3] + 2;
    cost_count[4] = cost_count[4] + 1;
    f = f + f.conjugate();
    f = Lucassequences(f,ceil((p+1)/fin));
    if flag == 1:
        f = 2;
    return f;

# Pairing (using new formulas and Jacobian coordinates)
# the order of the pairing is odd (the string has naf form)
def miller_Jacobian_odd_naf_trace_precomputed(X,Y,xx,yy,zz,AA,exp,fin):
    f=G(1);
    string = Trans(exp);
    l=len(string);
    flag=0;
    # Note that zz == 1
    zzs = 1;
    x1 = xx;
    y1 = yy;
    i = 1;
    while (i != l):
        [xx,yy,zz,tmp] = DBL_Jacobian(xx,yy,zz,zzs,AA);
        cost_count[2] = cost_count[2] - 2;
        cost_count[3] = cost_count[3] - 7;
        if string[i] == 1 and i != l-1:
            if xx != x1*zz^2:
                [xx2,yy2,zz2,Z1Z1,t0,H,r] = ADD_Jacobian(xx,yy,zz,x1,y1);
                cost_count[2] = cost_count[2] - 7;
                cost_count[3] = cost_count[3] - 4;
                f = line_fun_Jacobian_dbladd(xx,yy,zz,x1,y1,tmp,Z1Z1,t0,H,r,X,Y,f);
                xx = xx2;
                yy = yy2;
                zz = zz2;
                zzs = zz^2;
                #cost_count[3] = cost_count[3] + 1;
                i = i + 1;
            else:
                flag = 1;
        else:
            if string[i] == -1 and i != l-1:
                if xx != x1*zz^2:
                    [xx2,yy2,zz2,Z1Z1,t0,H,r] = ADD_Jacobian(xx,yy,zz,x1,-y1);
                    cost_count[2] = cost_count[2] - 7;
                    cost_count[3] = cost_count[3] - 4;
                    f = line_fun_Jacobian_dblmin(xx,yy,zz,x1,y1,tmp,Z1Z1,t0,H,r,X,Y,f);
                    xx = xx2;
                    yy = yy2;
                    zz = zz2;
                    zzs = zz^2;
                    #cost_count[3] = cost_count[3] + 1;
                    i = i + 1;
            else:
                if string[i] == 0 and string[i+1] == 0 and i < l-2:
                    zzs = zz^2;
                    #cost_count[3] = cost_count[3] + 1;
                    [xx2,yy2,zz2,tmp2] = DBL_Jacobian(xx,yy,zz,zzs,AA);
                    cost_count[2] = cost_count[2] - 2;
                    cost_count[3] = cost_count[3] - 7;
                    f = line_fun_Jacobian_qua(xx,yy,zz,zzs,X,Y,tmp,tmp2,f);
                    i = i + 2;
                    xx = xx2;
                    yy = yy2;
                    zz = zz2;
                    zzs = zz2^2;
                    #cost_count[3] = cost_count[3] + 1;
                else:
                    [f,zzs] = line_fun_double(xx,yy,zz,X,Y,tmp,f);
                    i = i + 1;
    if string[l-1] == 1:
        #if X != x1:
            #f = f * (X-x1);
            #cost_count[2] = cost_count[2] + 2;
        if X == x1:
            flag = 1;
    f = f^(p-1);
    cost_count[2] = cost_count[2] + 3;
    cost_count[3] = cost_count[3] + 2;
    cost_count[4] = cost_count[4] + 1;
    f = f + f.conjugate();
    f = Lucassequences(f,ceil((p+1)/fin));
    if flag == 1:
        f = 2;
    return f;

# Pairing (using new formulas and Jacobian coordinates)
# the order of the pairing is odd (the string has naf form)
def miller_Jacobian_odd_naf_trace_3(X,Y,xx,yy,zz,AA,exp,fin):
    f=G(1);
    string = Trans(exp);
    l=len(string);
    flag=0;
    # Note that zz == 1
    zzs = zz;
    x1 = xx;
    y1 = yy;
    [xx,yy,zz,tmp] = DBL_Jacobian(xx,yy,zz,zzs,AA);
    [f,zzs] = line_fun_double(xx,yy,zz,X,Y,tmp,f);
    f = f^(p-1);
    cost_count[2] = cost_count[2] + 3;
    cost_count[3] = cost_count[3] + 2;
    cost_count[4] = cost_count[4] + 1;
    f = f + f.conjugate();
    f = Lucassequences(f,ceil((p+1)/fin));
    if flag == 1:
        f = 2;
    return f;

# Pairing (using new formulas and Jacobian coordinates)
# the order of the pairing is odd (the string has naf form)
def miller_Jacobian_odd_naf_trace_3_precomputed(X,Y,xx,yy,zz,AA,exp,fin):
    f=G(1);
    string = Trans(exp);
    l=len(string);
    flag=0;
    # Note that zz == 1
    zzs = zz;
    x1 = xx;
    y1 = yy;
    [xx,yy,zz,tmp] = DBL_Jacobian(xx,yy,zz,zzs,AA);
    cost_count[2] = cost_count[2] - 2;
    cost_count[3] = cost_count[3] - 7;
    [f,zzs] = line_fun_double(xx,yy,zz,X,Y,tmp,f);
    f = f^(p-1);
    cost_count[2] = cost_count[2] + 3;
    cost_count[3] = cost_count[3] + 2;
    cost_count[4] = cost_count[4] + 1;
    f = f + f.conjugate();
    f = Lucassequences(f,ceil((p+1)/fin));
    if flag == 1:
        f = 2;
    return f;

# Batch cofactor exponentiation using trace
def batch_cof_trace(fi,f,ind,n, L):
	if n == 1:
		fi.append(f);
		return 1;
	m=floor(n/2);
	u=1;
	v=1;
	for i in range(m):
		u=L[ind[i]]*u;
	for j in range(m,n):
		v=L[ind[j]]*v;
	left=Lucassequences(f,u);
	right=Lucassequences(f,v);
	batch_cof_trace(fi,left,[ind[i] for i in range(m,n)],n-m, L);
	batch_cof_trace(fi,right,[ind[i] for i in range(m)],m, L);


AA = G(1-A^2/3); BB = G(A^3/27*2-A/3);
#y^2 = x^3+(1-A^2/9)x-A/3
EE = EllipticCurve(G,[AA,BB]);
inf = EE([0,1,0]);

Low = 10000000; High = 0;
cost_count = [0,0,0,0,0]; # M S m s i
print("Pairing computation in Jacobian coordinates:");
for count_time in range(10000):
    Plist = [];
    Qlist = [];
    tmp_cost_count = [];
    for i in range(5):
        tmp_cost_count.append(cost_count[i]);
    [P,Q]=elligator(A, E);
    [XP1, ZP1, XP2, ZP2] = Montgomery_ladder(4, A, P[0], 1);
    [XP,YP,ZP] = recover_y3(P[0], P[1], 1, XP1, ZP1, XP2, ZP2, A);
    [XQ1, ZQ1, XQ2, ZQ2] = Montgomery_ladder(4, A, Q[0], 1);
    [XQ,YQ,ZQ] = recover_y3(Q[0], Q[1], 1, XQ1, ZQ1, XQ2, ZQ2, A);
    ZPZQ = ZP*ZQ;
    ZPZQ = 1/ZPZQ;
    tZP = ZP;
    ZP = ZPZQ * ZQ;
    XP = XP * ZP;
    YP = YP * ZP;
    ZQ = ZPZQ * tZP;
    XQ = XQ * ZQ;
    YQ = YQ * ZQ;
    P = E([XP,YP]);
    Q = E([XQ,YQ]);
    cost_count[2] = cost_count[2] + 7;
    cost_count[4] = cost_count[4] + 1;
    f = miller_Jacobian_odd_naf_trace(Q[0]+G(A/3),Q[1],XP+G(A/3),YP,1,AA,(p+1)/4,(p+1)/4);
    #f = miller_modified_odd_naf_trace(Q[0]+G(A/3),Q[1],XP+G(A/3),YP,1,AA,AA,(p+1)/4,(p+1)/4);
    ind = [];
    n = 74; L = ells;
    fi=[]; ind=[]; If =[];
    for ii in range(len(ells)):
        ind.append(ii);
    batch_cof_trace(fi,f,ind,len(ells),ells);
    Ip = []; Iq = [];
    for i in range(n-1,-1,-1):
        if fi[i] == 2:
            Ip.append(n-1-i); Iq.append(n-1-i);If.append(n-1-i);
    Ipq = Ip;
    prime_prod = 1;
    for i in list(set(ind)-set(Ip)):
        prime_prod = prime_prod * ells[i];
    PXi = []; PZi = [];
    XPP = 0; YPP = 1; ZPP = 0;
    PP = E([0,1,0]);
    if Ip != []:
        R0X, R0Z, R1X, R1Z = Montgomery_ladder(prime_prod,A,P[0], 1);
        batch_cof_mul(PXi,PZi,R0X,R0Z,Ip,len(Ip),A,ells);
        if R0Z != 0:
            [XPP, YPP, ZPP] = recover_y3(P[0], P[1], 1, R0X, R0Z, R1X, R1Z, A);
            XPP = XPP + G(A/3)*ZPP;
            cost_count[2] = cost_count[2] + 1;
        else:
            XPP = 0; YPP = 1; ZPP = 0;
    tmp=Ip; Ip = [];
    for i in range(len(tmp)):
        if PZi[len(tmp)-1-i] == 0:
            Ip.append(tmp[i]);
    R0X, R0Z, R1X, R1Z = Montgomery_ladder(prime_prod,A,Q[0], 1);
    Q1Xi = []; Q1Zi = [];
    if Iq != []:
        batch_cof_mul(Q1Xi,Q1Zi,R0X,R0Z,Iq,len(Iq),A,ells);
    tmp=Iq; Iq = [];
    for i in range(len(tmp)):
        if Q1Zi[len(tmp)-1-i] == 0:
            Iq.append(tmp[i]);
    XP = P[0]; YP = P[1];
    XP = XP + G(A/3);
    ZP = 1;
    P = EE([XP,YP]);
    while Ip != []:
        [P1,Q1]=elligator(A, E);
        Qlist.append(Q1);
        prime_prod = 1;
        for i in list(set(ind)-set(Ip)):
            prime_prod = prime_prod * ells[i];
        R0X, R0Z, R1X, R1Z = Montgomery_ladder(4*prime_prod,A,P1[0], 1);
        [XP1,YP1,ZP1] = recover_y3(P1[0],P1[1],1,R0X,R0Z,R1X,R1Z,A);
        if ZP1 != 0:
            #ZP1 = 1/ZP1;
            #XP1 = XP1 * ZP1;
            #YP1 = YP1 * ZP1;
            #P1 = E([XP1,YP1]);
            P1Xi = []; P1Zi = [];
            batch_cof_mul(P1Xi,P1Zi,XP1,ZP1,Ip,len(Ip),A,ells);
            tmp=Ip; Ip = [];  scalar = 1;
            for i in range(len(tmp)):
                if P1Zi[len(tmp)-1-i] == 0:
                    Ip.append(tmp[i]);
                    scalar = scalar * ells[tmp[i]];
                else:
                    for ii in range(len(Ipq)):
                        if tmp[i] == Ipq[ii]:
                            PXi[len(Ipq)-1-ii] = P1Xi[len(tmp)-1-i];
                            PZi[len(Ipq)-1-ii] = P1Zi[len(tmp)-1-i];
            #R0X, R0Z, R1X, R1Z = Montgomery_ladder(scalar,A,XP1,ZP1);
            #[XP1,YP1,ZP1] = recover_y3(XP1,YP1,ZP1,R0X,R0Z,R1X,R1Z,A);
            if ZP1 != 0:
                #ZP1 = 1/ZP1;
                #XP1 = XP1 * ZP1;
                #YP1 = YP1 * ZP1;
                #XP1 = XP1 + G(A/3);
                #P1 = EE([XP1,YP1]);
                #[XP,YP,ZP] = ADD_Projective(XP,YP,ZP,XP1,YP1,1,AA,BB);
                [XP,YP,ZP] = ADD_Projective(XP,YP,ZP,XP1+G(A/3)*ZP1,YP1,ZP1,AA,BB);
                cost_count[2] = cost_count[2] + 1;
                if ZPP == 0:
                    XPP = XP1 + G(A/3)*ZP1; YPP = YP1; ZPP = ZP1;
                else:
                    [XPP, YPP, ZPP] = ADD_Projective(XPP,YPP,ZPP,XP1+G(A/3)*ZP1,YP1,ZP1,AA,BB);
            else:
                Ip = tmp;
    if ZP != 1:
        ZP = 1/G(ZP);
        XP = XP*ZP;
        YP = YP*ZP;
        cost_count[2] = cost_count[2] + 2;
        cost_count[4] = cost_count[4] + 1;
    XP = XP - G(A/3);
    P = E([XP,YP]);
    #PXi = []; PZi = []; ind = [];
    #n = 74;
    #for i in range(n):
    #    ind.append(i);
    #batch_cof_mul(PXi,PZi,XP,1,ind,n,A,ells);
    #Ip=[];
    #for i in range(n):
    #    if PZi[i] == 0:
    #        Ip.append(n-1-i);
    #print("a",Ip);
    #for i in range(n):
    #    xx,zz,xxx,zzz = Montgomery_ladder(ells[n-1-i],PXi[i],PZi[i],A);
    #    if zz == 0:
    #        print(i);
    XP = XP + G(A/3);
    [XQ,YQ] = [Q[0],Q[1]];
    ZQ = 1; YQ = Q[1];
    XQ = XQ + G(A/3);
    Q = EE([XQ,YQ]);
    try_time = 1;
    scalar = 1;
    if ZPP != 0:
        #ZPP = 1/ZPP;
        #XPP = XPP * ZPP;
        #YPP = YPP * ZPP;
        XPP = XPP - G(A/3)*ZPP;
        cost_count[2] = cost_count[2] + 1;
        #PP = E([XPP,YPP,ZPP]);
        scalar = 1;
        if Iq != If:
            if Iq != []:
                for i in list(set(If)-set(Iq)):
                    scalar = scalar * ells[i];
                R0X,R0Z,R1X,R1Z = Montgomery_ladder(scalar,A,XPP,ZPP);
                [XPP,YPP,ZPP] = recover_y3(XPP,YPP,ZPP,R0X,R0Z,R1X,R1Z,A);
        XPP = XPP + G(A/3)*ZPP;
        #PP = EE([XPP,YPP,ZPP]);
        ZPP = 1/ZPP;
        XPP = XPP * ZPP;
        YPP = YPP * ZPP;
        cost_count[2] = cost_count[2] + 2;
        cost_count[4] = cost_count[4] + 1;
        PP = EE([XPP,YPP]);
    else:
        if Iq != []:
            print("PP is error");
    while Iq != []:
        if try_time <= len(Qlist):
            Q1 = Qlist[try_time-1]; try_time = try_time + 1;
        else:
            [P1,Q1]=elligator(A, E);
        XQ1 = Q1[0]; YQ1 = Q1[1];
        prime_prod = 1;
        for i in list(set(ind)-set(Iq)):
            prime_prod = prime_prod * ells[i];
        remain = ceil((p+1)/(4*prime_prod));
        #print("tmp",PP.order(),remain,Iq,If);
        XQ1 = XQ1 + G(A/3); ZQ1 = 0;
        tmp=Iq; Iq = []; scalar = 1;
        if(tmp == [0]):
            if try_time == 2:
                f1 = miller_Jacobian_odd_naf_trace_3(XQ1,YQ1,XPP,YPP,1,AA,remain,remain);
            else:
                f1 = miller_Jacobian_odd_naf_trace_3_precomputed(XQ1,YQ1,XPP,YPP,1,AA,remain,remain);
        else:
            if try_time == 2:
                f1 = miller_Jacobian_odd_naf_trace(XQ1,YQ1,XPP,YPP,1,AA,remain,remain);
            else:
                f1 = miller_Jacobian_odd_naf_trace_precomputed(XQ1,YQ1,XPP,YPP,1,AA,remain,remain);
        f1i = []; If1 = 0;
        batch_cof_trace(f1i,f1,tmp,len(tmp),ells);
        for i in range(len(f1i)-1,-1,-1):
            if f1i[i] == 2:
                If1 = 1;
        if If1 == 0:
            XQ1 = XQ1 - G(A/3);
            R0X, R0Z, R1X, R1Z = Montgomery_ladder(4*prime_prod,A,XQ1,1);
            [XQ1,YQ1,ZQ1] = recover_y3(Q1[0],Q1[1],1,R0X,R0Z,R1X,R1Z,A);
        if ZQ1 != 0:
            #ZQ1 = 1/ZQ1;
            #XQ1 = XQ1 * ZQ1;
            #YQ1 = YQ1 * ZQ1;
            #XQ1 = XQ1 + G(A/3);
            #Q1 = EE([XQ1,YQ1]);
            #[XQ,YQ,ZQ] = ADD_Projective(XQ,YQ,ZQ,XQ1+G(A/3),YQ1,1,AA,BB);
            [XQ,YQ,ZQ] = ADD_Projective(XQ,YQ,ZQ,XQ1+G(A/3)*ZQ1,YQ1,ZQ1,AA,BB);
            cost_count[2] = cost_count[2] + 1;
        else:
            Iq = tmp;
    #Q = EE([XQ/ZQ,YQ/ZQ]);
    if ZQ != 1:
        ZQ = 1/G(ZQ);
        XQ = XQ * ZQ;
        YQ = YQ * ZQ;
        cost_count[2] = cost_count[2] + 2;
        cost_count[4] = cost_count[4] + 1;
    XQ = XQ-G(A/3);
    Q = E([XQ,YQ]);
    for i in range(5):
        tmp_cost_count[i] = cost_count[i] - tmp_cost_count[i];
    cost_sum = tmp_cost_count[0]*3+tmp_cost_count[1]*2+tmp_cost_count[2]*1+tmp_cost_count[3]*0.8+tmp_cost_count[4]*30;
    if Low >= cost_sum:
        Low = cost_sum;
    if High <= cost_sum:
        High = cost_sum;
    #QXi = []; QZi = []; ind = [];
    #n = 74;
    #for i in range(n):
    #    ind.append(i);
    #batch_cof_mul(QXi,QZi,XQ,1,ind,n,A,ells);
    #Iq=[];
    #for i in range(n):
    #    if QZi[i] == 0:
    #        Iq.append(n-1-i);
    #print(count_time,Iq);
    #for i in range(n):
    #    xx,zz,xxx,zzz = Montgomery_ladder(ells[n-1-i],QXi[i],QZi[i],A);
    #    if zz == 0:
    #        print(i);

#if Ip == [] and Iq == []:
#    print("Torsion basis are generated");
#else:
#    print("error");
#print(cost_count);
#print(cost_count[0]*3+cost_count[1]*2+cost_count[2]*1+cost_count[3]*0.8+cost_count[4]*30)
print("Average cost:",(cost_count[0]*3+cost_count[1]*2+cost_count[2]*1+cost_count[3]*0.8+cost_count[4]*30)/(count_time+1));
print("Lowest:",Low,"Highest:",High);

Low = 10000000; High = 0;
cost_count = [0,0,0,0,0]; # M S m s i
print("Pairing computation in modified Jacobian coordinates:");
bench = 10000
print("Bench = ",bench);
for count_time in range(bench):
    Plist = [];
    Qlist = [];
    tmp_cost_count = [];
    for i in range(5):
        tmp_cost_count.append(cost_count[i]);
    [P,Q]=elligator(A, E);
    [XP1, ZP1, XP2, ZP2] = Montgomery_ladder(4, A, P[0], 1);
    [XP,YP,ZP] = recover_y3(P[0], P[1], 1, XP1, ZP1, XP2, ZP2, A);
    [XQ1, ZQ1, XQ2, ZQ2] = Montgomery_ladder(4, A, Q[0], 1);
    [XQ,YQ,ZQ] = recover_y3(Q[0], Q[1], 1, XQ1, ZQ1, XQ2, ZQ2, A);
    ZPZQ = ZP*ZQ;
    ZPZQ = 1/ZPZQ;
    tZP = ZP;
    ZP = ZPZQ * ZQ;
    XP = XP * ZP;
    YP = YP * ZP;
    ZQ = ZPZQ * tZP;
    XQ = XQ * ZQ;
    YQ = YQ * ZQ;
    P = E([XP,YP]);
    Q = E([XQ,YQ]);
    cost_count[2] = cost_count[2] + 7;
    cost_count[4] = cost_count[4] + 1;
    #f = miller_Jacobian_odd_naf_trace(Q[0]+G(A/3),Q[1],XP+G(A/3),YP,1,AA,(p+1)/4,(p+1)/4);
    f = miller_modified_odd_naf_trace(Q[0]+G(A/3),Q[1],XP+G(A/3),YP,1,AA,AA,(p+1)/4,(p+1)/4);
    ind = [];
    n = 74; L = ells;
    fi=[]; ind=[]; If =[];
    for ii in range(len(ells)):
        ind.append(ii);
    batch_cof_trace(fi,f,ind,len(ells),ells);
    Ip = []; Iq = [];
    for i in range(n-1,-1,-1):
        if fi[i] == 2:
            Ip.append(n-1-i); Iq.append(n-1-i);If.append(n-1-i);
    Ipq = Ip;
    prime_prod = 1;
    for i in list(set(ind)-set(Ip)):
        prime_prod = prime_prod * ells[i];
    PXi = []; PZi = [];
    XPP = 0; YPP = 1; ZPP = 0;
    PP = E([0,1,0]);
    if Ip != []:
        R0X, R0Z, R1X, R1Z = Montgomery_ladder(prime_prod,A,P[0], 1);
        batch_cof_mul(PXi,PZi,R0X,R0Z,Ip,len(Ip),A,ells);
        if R0Z != 0:
            [XPP, YPP, ZPP] = recover_y3(P[0], P[1], 1, R0X, R0Z, R1X, R1Z, A);
            XPP = XPP + G(A/3)*ZPP;
            cost_count[2] = cost_count[2] + 1;
        else:
            XPP = 0; YPP = 1; ZPP = 0;
    tmp=Ip; Ip = [];
    for i in range(len(tmp)):
        if PZi[len(tmp)-1-i] == 0:
            Ip.append(tmp[i]);
    R0X, R0Z, R1X, R1Z = Montgomery_ladder(prime_prod,A,Q[0], 1);
    Q1Xi = []; Q1Zi = [];
    if Iq != []:
        batch_cof_mul(Q1Xi,Q1Zi,R0X,R0Z,Iq,len(Iq),A,ells);
    tmp=Iq; Iq = [];
    for i in range(len(tmp)):
        if Q1Zi[len(tmp)-1-i] == 0:
            Iq.append(tmp[i]);
    XP = P[0]; YP = P[1];
    XP = XP + G(A/3);
    ZP = 1;
    P = EE([XP,YP]);
    while Ip != []:
        [P1,Q1]=elligator(A, E);
        Qlist.append(Q1);
        prime_prod = 1;
        for i in list(set(ind)-set(Ip)):
            prime_prod = prime_prod * ells[i];
        R0X, R0Z, R1X, R1Z = Montgomery_ladder(4*prime_prod,A,P1[0], 1);
        [XP1,YP1,ZP1] = recover_y3(P1[0],P1[1],1,R0X,R0Z,R1X,R1Z,A);
        if ZP1 != 0:
            #ZP1 = 1/ZP1;
            #XP1 = XP1 * ZP1;
            #YP1 = YP1 * ZP1;
            #P1 = E([XP1,YP1]);
            P1Xi = []; P1Zi = [];
            batch_cof_mul(P1Xi,P1Zi,XP1,ZP1,Ip,len(Ip),A,ells);
            tmp=Ip; Ip = [];  scalar = 1;
            for i in range(len(tmp)):
                if P1Zi[len(tmp)-1-i] == 0:
                    Ip.append(tmp[i]);
                    scalar = scalar * ells[tmp[i]];
                else:
                    for ii in range(len(Ipq)):
                        if tmp[i] == Ipq[ii]:
                            PXi[len(Ipq)-1-ii] = P1Xi[len(tmp)-1-i];
                            PZi[len(Ipq)-1-ii] = P1Zi[len(tmp)-1-i];
            #R0X, R0Z, R1X, R1Z = Montgomery_ladder(scalar,A,XP1,ZP1);
            #[XP1,YP1,ZP1] = recover_y3(XP1,YP1,ZP1,R0X,R0Z,R1X,R1Z,A);
            if ZP1 != 0:
                #ZP1 = 1/ZP1;
                #XP1 = XP1 * ZP1;
                #YP1 = YP1 * ZP1;
                #XP1 = XP1 + G(A/3);
                #P1 = EE([XP1,YP1]);
                #[XP,YP,ZP] = ADD_Projective(XP,YP,ZP,XP1,YP1,1,AA,BB);
                [XP,YP,ZP] = ADD_Projective(XP,YP,ZP,XP1+G(A/3)*ZP1,YP1,ZP1,AA,BB);
                cost_count[2] = cost_count[2] + 1;
                if ZPP == 0:
                    XPP = XP1 + G(A/3)*ZP1; YPP = YP1; ZPP = ZP1;
                else:
                    [XPP, YPP, ZPP] = ADD_Projective(XPP,YPP,ZPP,XP1+G(A/3)*ZP1,YP1,ZP1,AA,BB);
            else:
                Ip = tmp;
    if ZP != 1:
        ZP = 1/G(ZP);
        XP = XP*ZP;
        YP = YP*ZP;
        cost_count[2] = cost_count[2] + 2;
        cost_count[4] = cost_count[4] + 1;
    XP = XP - G(A/3);
    P = E([XP,YP]);
    #PXi = []; PZi = []; ind = [];
    #n = 74;
    #for i in range(n):
    #    ind.append(i);
    #batch_cof_mul(PXi,PZi,XP,1,ind,n,A,ells);
    #Ip=[];
    #for i in range(n):
    #    if PZi[i] == 0:
    #        Ip.append(n-1-i);
    #print("a",Ip);
    #for i in range(n):
    #    xx,zz,xxx,zzz = Montgomery_ladder(ells[n-1-i],PXi[i],PZi[i],A);
    #    if zz == 0:
    #        print(i);
    XP = XP + G(A/3);
    [XQ,YQ] = [Q[0],Q[1]];
    ZQ = 1; YQ = Q[1];
    XQ = XQ + G(A/3);
    Q = EE([XQ,YQ]);
    try_time = 1;
    scalar = 1;
    if ZPP != 0:
        #ZPP = 1/ZPP;
        #XPP = XPP * ZPP;
        #YPP = YPP * ZPP;
        XPP = XPP - G(A/3)*ZPP;
        cost_count[2] = cost_count[2] + 1;
        #PP = E([XPP,YPP,ZPP]);
        scalar = 1;
        if Iq != If:
            if Iq != []:
                for i in list(set(If)-set(Iq)):
                    scalar = scalar * ells[i];
                R0X,R0Z,R1X,R1Z = Montgomery_ladder(scalar,A,XPP,ZPP);
                [XPP,YPP,ZPP] = recover_y3(XPP,YPP,ZPP,R0X,R0Z,R1X,R1Z,A);
        XPP = XPP + G(A/3)*ZPP;
        #PP = EE([XPP,YPP,ZPP]);
        ZPP = 1/ZPP;
        XPP = XPP * ZPP;
        YPP = YPP * ZPP;
        cost_count[2] = cost_count[2] + 2;
        cost_count[4] = cost_count[4] + 1;
        PP = EE([XPP,YPP]);
    else:
        if Iq != []:
            print("PP is error");
    while Iq != []:
        if try_time <= len(Qlist):
            Q1 = Qlist[try_time-1]; try_time = try_time + 1;
        else:
            [P1,Q1]=elligator(A, E);
        XQ1 = Q1[0]; YQ1 = Q1[1];
        prime_prod = 1;
        for i in list(set(ind)-set(Iq)):
            prime_prod = prime_prod * ells[i];
        remain = ceil((p+1)/(4*prime_prod));
        #print("tmp",PP.order(),remain,Iq,If);
        XQ1 = XQ1 + G(A/3); ZQ1 = 0;
        tmp=Iq; Iq = []; scalar = 1;
        if(tmp == [0]):
            if try_time == 2:
                f1 = miller_Jacobian_odd_naf_trace_3(XQ1,YQ1,XPP,YPP,1,AA,remain,remain);
            else:
                f1 = miller_Jacobian_odd_naf_trace_3_precomputed(XQ1,YQ1,XPP,YPP,1,AA,remain,remain);
        else:
            if try_time == 2:
                f1 = miller_Jacobian_odd_naf_trace(XQ1,YQ1,XPP,YPP,1,AA,remain,remain);
            else:
                f1 = miller_Jacobian_odd_naf_trace_precomputed(XQ1,YQ1,XPP,YPP,1,AA,remain,remain);
        f1i = []; If1 = 0;
        batch_cof_trace(f1i,f1,tmp,len(tmp),ells);
        for i in range(len(f1i)-1,-1,-1):
            if f1i[i] == 2:
                If1 = 1;
        if If1 == 0:
            XQ1 = XQ1 - G(A/3);
            R0X, R0Z, R1X, R1Z = Montgomery_ladder(4*prime_prod,A,XQ1,1);
            [XQ1,YQ1,ZQ1] = recover_y3(Q1[0],Q1[1],1,R0X,R0Z,R1X,R1Z,A);
        if ZQ1 != 0:
            #ZQ1 = 1/ZQ1;
            #XQ1 = XQ1 * ZQ1;
            #YQ1 = YQ1 * ZQ1;
            #XQ1 = XQ1 + G(A/3);
            #Q1 = EE([XQ1,YQ1]);
            #[XQ,YQ,ZQ] = ADD_Projective(XQ,YQ,ZQ,XQ1+G(A/3),YQ1,1,AA,BB);
            [XQ,YQ,ZQ] = ADD_Projective(XQ,YQ,ZQ,XQ1+G(A/3)*ZQ1,YQ1,ZQ1,AA,BB);
            cost_count[2] = cost_count[2] + 1;
        else:
            Iq = tmp;
    #Q = EE([XQ/ZQ,YQ/ZQ]);
    if ZQ != 1:
        ZQ = 1/G(ZQ);
        XQ = XQ * ZQ;
        YQ = YQ * ZQ;
        cost_count[2] = cost_count[2] + 2;
        cost_count[4] = cost_count[4] + 1;
    XQ = XQ-G(A/3);
    Q = E([XQ,YQ]);
    for i in range(5):
        tmp_cost_count[i] = cost_count[i] - tmp_cost_count[i];
    cost_sum = tmp_cost_count[0]*3+tmp_cost_count[1]*2+tmp_cost_count[2]*1+tmp_cost_count[3]*0.8+tmp_cost_count[4]*30;
    if Low >= cost_sum:
        Low = cost_sum;
    if High <= cost_sum:
        High = cost_sum;
    #QXi = []; QZi = []; ind = [];
    #n = 74;
    #for i in range(n):
    #    ind.append(i);
    #batch_cof_mul(QXi,QZi,XQ,1,ind,n,A,ells);
    #Iq=[];
    #for i in range(n):
    #    if QZi[i] == 0:
    #        Iq.append(n-1-i);
    #print(count_time,Iq);
    #for i in range(n):
    #    xx,zz,xxx,zzz = Montgomery_ladder(ells[n-1-i],QXi[i],QZi[i],A);
    #    if zz == 0:
    #        print(i);

#if Ip == [] and Iq == []:
#    print("Torsion basis are generated");
#else:
#    print("error");
#print(cost_count);
#print(cost_count[0]*3+cost_count[1]*2+cost_count[2]*1+cost_count[3]*0.8+cost_count[4]*30)
print("Average cost:",(cost_count[0]*3+cost_count[1]*2+cost_count[2]*1+cost_count[3]*0.8+cost_count[4]*30)/(count_time+1));
print("Lowest",Low,"Highest",High);