
test_case =6;
Nx=300;
Ny=300;

output_name = fullfile('outputs',[output,num2str(test_case)]);
[Data] = spec2meshgrid_flat_shell(test_case,Nx,Ny,'displacement',3,'upper',output_name); % Uz 

[Data] = spec2meshgrid_flat_shell(test_case,Nx,Ny,'displacement',8,'upper',output_name); % sqrt((Ux+h/2.*Fix).^2+(Uy+h/2.*Fiy).^2)