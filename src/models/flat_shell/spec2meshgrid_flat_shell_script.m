
test_case =6;
Nx=300;
Ny=300;

[Data] = spec2meshgrid_flat_shell(test_case,Nx,Ny,'displacement',3,'upper'); % Uz 

[Data] = spec2meshgrid_flat_shell(test_case,Nx,Ny,'displacement',8,'upper'); % sqrt((Ux+h/2.*Fix).^2+(Uy+h/2.*Fiy).^2)