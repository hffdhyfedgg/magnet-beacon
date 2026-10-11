function rm = CalcRotMat(ax, ay, az, bx, by, bz, u, v, w, r, b)



cx = 3*b.*u.*v./r.^3;
cy = b.*(3*v.^2 - 1)./r.^3;
cz = 3*b.*v.*w./r.^3;
dx = b.*(3*u.^2 - 1)./r.^3;
dy = 3*b.*u.*v./r.^3;
dz = 3*b.*u.*w./r.^3;

var = ((cy.^2 + cz.^2).*dx.^2 + (- 2.*cx.*cy.*dy - 2.*cx.*cz.*dz).*dx + (cx.^2.*dy.^2 + cx.^2.*dz.^2 + cy.^2.*dz.^2 - 2.*cy.*cz.*dy.*dz + cz.^2.*dy.^2));



r_11 = ((bx.*cy.^2 - ax.*dy.*cy + bx.*cz.^2 - ax.*dz.*cz).*dx + (ax.*cx.*dy.^2 + ax.*cx.*dz.^2 + ay.*bz.*cy.*dz - ay.*bz.*cz.*dy - az.*by.*cy.*dz + az.*by.*cz.*dy - bx.*cx.*cy.*dy - bx.*cx.*cz.*dz))./  var;
r_12 = ((ax.*cy).*dx.^2 + (ay.*bz.*cz - az.*by.*cz - ax.*cx.*dy - bx.*cx.*cy).*dx + (ax.*cy.*dz.^2 + bx.*cx.^2.*dy + bx.*cz.^2.*dy - ay.*bz.*cx.*dz + az.*by.*cx.*dz - ax.*cz.*dy.*dz - bx.*cy.*cz.*dz))./var;
r_13 = ((ax.*cz).*dx.^2 + (az.*by.*cy - ay.*bz.*cy - ax.*cx.*dz - bx.*cx.*cz).*dx + (ax.*cz.*dy.^2 + bx.*cx.^2.*dz + bx.*cy.^2.*dz + ay.*bz.*cx.*dy - az.*by.*cx.*dy - ax.*cy.*dy.*dz - bx.*cy.*cz.*dy))./var;
r_21 = ((by.*cy.^2 - ay.*dy.*cy + by.*cz.^2 - ay.*dz.*cz).*dx + (ay.*cx.*dy.^2 + ay.*cx.*dz.^2 - ax.*bz.*cy.*dz + ax.*bz.*cz.*dy + az.*bx.*cy.*dz - az.*bx.*cz.*dy - by.*cx.*cy.*dy - by.*cx.*cz.*dz))./  var;
r_22 = ((ay.*cy).*dx.^2 + (az.*bx.*cz - ax.*bz.*cz - ay.*cx.*dy - by.*cx.*cy).*dx + (ay.*cy.*dz.^2 + by.*cx.^2.*dy + by.*cz.^2.*dy + ax.*bz.*cx.*dz - az.*bx.*cx.*dz - ay.*cz.*dy.*dz - by.*cy.*cz.*dz))./var;
r_23 = ((ay.*cz).*dx.^2 + (ax.*bz.*cy - az.*bx.*cy - ay.*cx.*dz - by.*cx.*cz).*dx + (ay.*cz.*dy.^2 + by.*cx.^2.*dz + by.*cy.^2.*dz - ax.*bz.*cx.*dy + az.*bx.*cx.*dy - ay.*cy.*dy.*dz - by.*cy.*cz.*dy))./var;
r_31 = ((bz.*cy.^2 - az.*dy.*cy + bz.*cz.^2 - az.*dz.*cz).*dx + (az.*cx.*dy.^2 + az.*cx.*dz.^2 + ax.*by.*cy.*dz - ax.*by.*cz.*dy - ay.*bx.*cy.*dz + ay.*bx.*cz.*dy - bz.*cx.*cy.*dy - bz.*cx.*cz.*dz))./  var;
r_32 = ((az.*cy).*dx.^2 + (ax.*by.*cz - ay.*bx.*cz - az.*cx.*dy - bz.*cx.*cy).*dx + (az.*cy.*dz.^2 + bz.*cx.^2.*dy + bz.*cz.^2.*dy - ax.*by.*cx.*dz + ay.*bx.*cx.*dz - az.*cz.*dy.*dz - bz.*cy.*cz.*dz))./var;
r_33 = ((az.*cz).*dx.^2 + (ay.*bx.*cy - ax.*by.*cy - az.*cx.*dz - bz.*cx.*cz).*dx + (az.*cz.*dy.^2 + bz.*cx.^2.*dz + bz.*cy.^2.*dz + ax.*by.*cx.*dy - ay.*bx.*cx.*dy - az.*cy.*dy.*dz - bz.*cy.*cz.*dy))./var;

rm = [r_11,r_12,r_13; r_21,r_22,r_23; r_31,r_32,r_33];
end

