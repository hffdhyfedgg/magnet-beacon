
b = 1;
x = 1;
y = 1;
z = 1;
r = sqrt(sum([x,y,z].^2));
u = x/r;
v = y/r;
w = z/r;
alf = deg2rad(23);
bet = deg2rad(34);
gam = deg2rad(56);

rot_mat = Ang2Mat(alf,bet,gam);
b_r = GenV_o(b, alf,bet, gam, y,x,z);



rm = CalcRotMat(repmat(b_r(1),1,1,3), repmat(b_r(2),1,1,3), repmat(b_r(3),1,1,3), repmat(b_r(4),1,1,3), repmat(b_r(5),1,1,3), repmat(b_r(6),1,1,3), repmat(u,1,1,3), repmat(v,1,1,3), repmat(w,1,1,3), repmat(r,1,1,3), repmat(b,1,1,3) );

angles = rad2deg(rotm2eul(rm, 'ZYX'))
