syms a b c d a11 a12 a13 a14 a21 a22 a23 a24

aa = [a,0,0,b;
     0,0,0,0;
     0,0,0,0;
     c,0,0,d];

bb= aa*aa.';
bb

aa = [a11,a12,a13,0;
      0,1,0,0;
      a21,0,a23,a24;
      0,0,0,1];
  
bb = kron(aa,aa.');
bb
size(bb)