function K = mod_Besselfunktion(u)
    z = (3*u/2).^(2/3);
    K = pi .* sqrt(3./z) .* airy(z) - (10*pi./(3.*u)) .* sqrt(3./z) .* airy(1,z);

     
end


x = 0:100:100000;
ai = airy(x);
bi = airy(1,x);
K = mod_Besselfunktion(x)
figure
plot(x,K,'-b',x,bi,'-r')
axis([0 100000 0 10^17])
xlabel('x')
legend('K5/3(x)','d/dx Ai(x)','Location','NorthWest')