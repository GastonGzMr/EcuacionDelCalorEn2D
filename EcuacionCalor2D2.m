clear; close all; clc;

%% Input
                                  
CONDUCTIVIDAD = 204.3;                                                     % Conductividad (J/m*s)                                     
CALOR_ESPECIFICO = 910;                                                    % Calor específico (J/Kg*K)
DENSIDAD = 2700.0;                                                         % Densidad (Kg/m^3)
hc = 200;                                                                  % Transmitancia térmica (W/m^2*K)                                                                  

L= 0.06;                                                                   % Longitud (m)
H= 0.04;                                                                   % Altura (m)
dx=0.005;                                                                  % Delta x (m)

t_fin= 20;                                                                 % Tiempo de evaluación (s)
dt= 0.05;                                                                  % Delta t (s)
  
T_ambiente= 293.15 ;                                                       % Temperatura mínima (K)
T_superficie = 343.15;                                                     % Temperatura máxima (K)
f=@(t) T_ambiente+(T_superficie-T_ambiente)*sin(t/10*pi);                  % Función de temperatura en el tiempo (K)


%% Constantes

Nx=L/dx;                                                                   % Nodos X
Ny=H/dx;                                                                   % Nodos Y
n_time=round(t_fin/dt);                                                    % Numero de iteraciones por tiempo
c= (hc*dx)/CONDUCTIVIDAD;                                                  % Coeficiente de la salida de calor por aire
r=(dt*CONDUCTIVIDAD)/(DENSIDAD*CALOR_ESPECIFICO*(dx^2));                   % Coeficiente de la transferencia de calor entre nodos

if r<= 1/4                                                                 % Verifica la condición de estabilidad.
else
    fprintf('Error, la condición de estabilidad no se cumple');
    return
end

%% Condiciones iniciales

M=zeros(Nx*Ny);
Z=zeros(2*Nx+1,1);                                                         % Defino el vector Z para construir la matriz M con mayor facilidad.
Z(1)=r;
Z(Nx)=r;
Z(1+Nx)=(1-4*r);
Z(2+Nx)=r;
Z(2*Nx+1)=r;
U=zeros(Nx*(Ny+1)-1,n_time);

%% Construcción de M
for i=(Nx+1):Nx*Ny
    if(mod(i-1,Nx)>0)
        if(mod(i,Nx)==0)
            M(i,:)=M(i-1,:);
        else
            M(i,(i-Nx):(i+Nx))=Z;                                          %Condición 2.
        end
    end
end

%% Condición 3.
M(2:Nx,:)= M(2+Nx:2*Nx,:); 

%% Condición 4.
Z(1+Nx)=(1-(3+c)*r);
Z(2*Nx+1)=c;
for i=Nx*Ny+1:Nx*(Ny-1)-1
    if(mod(i-1,Nx)>0)
        if(mod(i,Nx)==0)
            M(i,:)=M(i-1,:);
        else
            M(i,(i-Nx):(i+Nx))=Z;
        end
    end
end
F=zeros(Nx*Ny,n_time);
U(1:Nx*Ny,1)=T_ambiente;
U((Nx*Ny+1):Nx*(Ny+1)-1,:)=T_ambiente;

%% Condición 1
for i=1:n_time
    if i<=200
        F(1:Nx:Nx*Ny,i)=feval(f,i*dt);
    else
        F(1:Nx:Nx*Ny,i)=T_ambiente;
    end
end

%% Cálculo de U
for i=2:n_time
    U(1:Nx*Ny,i)=M*U(:,i-1)+F(:,i);
end

%% Plot


T=zeros(Ny,Nx,n_time);
for k=1:n_time
    for i=1:Ny
        T(i,1:Nx,k)=U((i-1)*Nx+1:i*Nx,k);
    end
end
        

x=zeros(1,Nx-1);y=zeros(1,Ny);                                             %Genera la superficie
for i = 1:Nx-1         
    x(i) =(i-1)*dx; 
end
for i = 1:Ny            
    y(i) =(i-1)*dx; 
end

eda= T(:,1:Nx-1,1);
xlim([0 L+dx]); ylabel('Largo');
ylim([0 H+dx]); xlabel('Ancho');
zlim([0 T_superficie]); zlabel('Temperatura');
for j=1:n_time                                
    surf(x,y,T(:,1:Nx-1,j))
    %cb=colorbar;
    %caxis([T_ambiente 303.15]);
    view(0,90);
    title(sprintf('Temperatura en la iteración %i ',j))
    drawnow
end
