function niz2 = prosiri(niz1,k)

% niz1 je kolona!

pkg load communications;
pkg load signal;

temp = niz1.';
niz2 = repmat(temp,k,1);
niz2 = niz2(:); 
