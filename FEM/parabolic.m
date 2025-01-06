classdef parabolic < handle
    % classe per gestire un problema parabolico con laplaciano, senza sorgente, misto con condizioni di Dirichlet con dato della forma g(x,t)=zeta(t)*G(x) e con condizioni di Neumann omogenee
    
    properties
        coord % coordinate dei vertici
        topol % topologia degli elementi
        bound % vertici della porzione di Dirichlet del bordo con relativo valore di G
        dt % step temporale
        T % tempo finale
        zeta % funzione t->zeta(t) per le conizioni di Dirichlet
        dzeta % derivata di t->zeta(t)
        u0 % dato iniziale
        u % soluzione
    end

    methods
        function obj=parabolic(path,data)
            % costruttore della classe, inizializza le coordinate dei vertici, la topologia degli elementi e le condizioni di Dirichlet da file
            
            obj.coord=load(sprintf("%s/%s/%s.coord",path,data,data));
            obj.topol=load(sprintf("%s/%s/%s.topol",path,data,data));
            obj.bound=load(sprintf("%s/%s/%s.bound",path,data,data));
        end

        function k=computeKloc(obj,e)
            % calcola la stiffness matrix locale relativa all'elemento e

            X=[ones(3,1),obj.coord(obj.topol(e,:),:)]; % matrice X così come definita nel report (in realtà è la trasposta)
            A=det(X)*inv(X)'; % matrice dei cofattori di X da cui si possono determinare a, b e c
            Delta=0.5*det(X); % area dell'elemento
            b=A(:,2); % vettore dei b
            c=A(:,3); % vettore dei c
            k=sparse(1/(4*Delta)*(b*b'+c*c')); % risultato finale facendo il prodotto tensore b*b' e c*c'
        end

        function k=assemblyK(obj)
            % assembla la stiffness matrix

            k=sparse(size(obj.coord,1),size(obj.coord,1));

            for e=1:size(obj.topol,1)
                k(obj.topol(e,:),obj.topol(e,:))=k(obj.topol(e,:),obj.topol(e,:))+obj.computeKloc(e); %#ok<SPRIX>
            end
        end

        function m=computeMloc(obj,e)
            % calcola la mass matrix locale relativa all'elemento e

            C=[ones(3,1),obj.coord(obj.topol(e,:),:)]; % matrice C così come definita a lezione
            Delta=0.5*det(C); % area dell'elemento
            m=sparse(Delta/12*[2,1,1;1,2,1;1,1,2]); % calcolo del risultato finale
        end

        function m=assemblyM(obj)
            % assembla la mass matrix

            m=sparse(size(obj.coord,1),size(obj.coord,1));

            for e=1:size(obj.topol,1)
                m(obj.topol(e,:),obj.topol(e,:))=m(obj.topol(e,:),obj.topol(e,:))+obj.computeMloc(e); %#ok<SPRIX>
            end
        end

        function [B,C]=computeMatrices(obj,K,M)
            % calola le matrici B e C come definite nel report

            C=M/obj.dt-K/2;
            C(obj.bound(:,1),:)=0;
            C(:,obj.bound(:,1))=0;

            B=M/obj.dt+K/2;
            B(obj.bound(:,1),:)=0;
            B(:,obj.bound(:,1))=0;
            B(obj.bound(:,1),obj.bound(:,1))=speye(size(obj.bound,1));
        end

        function Phi=computePhi(obj,K,M,C,u,t)
            % calcola Phi come definito nel report

            Phi=C*u;
            Phi=Phi-K(:,obj.bound(:,1))*0.5*(obj.zeta(t)+obj.zeta(t+obj.dt))*obj.bound(:,2);
            Phi=Phi-M(:,obj.bound(:,1))*0.5*(obj.dzeta(t)+obj.dzeta(t+obj.dt))*obj.bound(:,2);
            Phi(obj.bound(:,1))=obj.zeta(t+obj.dt)*obj.bound(:,2);
        end

        function w=solve(obj,precond)
            % risoluzione del problema parabolico

            char=fprintf("Inizio la risoluzione del problema parabolico\n");

            char=char+fprintf("Sto assemblando le matrici\n");

            K=obj.assemblyK(); % matrice di stiffness
            M=obj.assemblyM(); % matrice di massa

            [B,C]=obj.computeMatrices(K,M); % matrici B ed C così come definite nel report
            
            % calcolo del precondizionatore di B
            char=char+fprintf("Sto precondizionando\n");
            if precond=="IC"
                L=ichol(B);
            elseif precond=="J"
                L=sparse(diag(diag(B)));
            else
                L=speye(size(B,2));
            end
            
            N=fix(obj.T/obj.dt); % numero di step temporali
            obj.u=zeros(size(B,2),N+1); % prealloco lo spazio per la soluzione

            obj.u(:,1)=obj.u0; % dato iniziale
            t=0;
            
            char=char+fprintf("Sto risolvendo l'ODE: %3.0f/100",0);

            for k=1:N
                fprintf("\b\b\b\b\b\b\b%3.0f/100",k/N*100);
                r=obj.computePhi(K,M,C,obj.u(:,k),t); % calcolo del membro di destra del sistema lineare
                [obj.u(:,k+1),~,~]=PCG(B,r,zeros(size(B,2),1),1e-8,1e4,L); % risoluzione del sistema lineare con PCG
                t=t+obj.dt; % avanzamento del tempo
            end

            w=obj.u; % ritorna la soluzione

            fprintf(repmat('\b',1,char));
        end

        function plot(obj,k,style)
        	% disegna la soluzione al tempo t=(k-1)*dt con stile per le linee definito da style

            clim([min(min(obj.u)),max(max(obj.u))])
        
            w=obj.u;
            w=w(:,k);
            patch("Faces",obj.topol,"Vertices",obj.coord,"FaceVertexCData",w,"FaceColor","interp","LineStyle",style);
            axis equal;
            colorbar;
            xlim([min(obj.coord(:,1)),max(obj.coord(:,1))]);
            ylim([min(obj.coord(:,2)),max(obj.coord(:,2))]);
            N=fix(obj.T/obj.dt);
            title(sprintf("t=%2.2f",(k-1)/N*obj.T));
        end

        function animate(obj,style)
        	% anima la soluzione con stile per le linee definito da style

            clim([min(min(obj.u)),max(max(obj.u))]);

            w=obj.u;
            patch("Faces",obj.topol,"Vertices",obj.coord,"LineStyle",style);
            axis equal;
            colorbar;
            xlim([min(obj.coord(:,1)),max(obj.coord(:,1))]);
            ylim([min(obj.coord(:,2)),max(obj.coord(:,2))]);
            hold on;
            p=patch("Faces",obj.topol,"Vertices",obj.coord,"FaceVertexCData",w(:,1),"FaceColor","interp","LineStyle",style);
            hold off;

            N=fix(obj.T/obj.dt);

            for k=1:N+1
                p.FaceVertexCData=w(:,k);
                title(sprintf("t=%2.2f",(k-1)/N*obj.T));
                drawnow;
            end
        end

        function movie(obj,style,filename)
        	% crea video dell'animazione della soluzione con stile per le linee definito da style

            clim([min(min(obj.u)),max(max(obj.u))])

            w=obj.u;
            patch("Faces",obj.topol,"Vertices",obj.coord,"LineStyle",style);
            axis equal;
            colorbar;
            xlim([min(obj.coord(:,1)),max(obj.coord(:,1))]);
            ylim([min(obj.coord(:,2)),max(obj.coord(:,2))]);
            hold on;
            p=patch("Faces",obj.topol,"Vertices",obj.coord,"FaceVertexCData",w(:,1),"FaceColor","interp","LineStyle",style);
            hold off;

            N=fix(obj.T/obj.dt);

            v=VideoWriter(sprintf("%s.avi",filename));
            v.FrameRate=7;
            v.Quality=100;
            open(v);

            for k=1:N+1
                fprintf("%3.0f/100",k/(N+1)*100);
                p.FaceVertexCData=w(:,k);
                title(sprintf("t=%2.2f",(k-1)/N*obj.T));
                writeVideo(v,getframe);
                fprintf("\b\b\b\b\b\b\b");
            end

            fprintf("\n");

            close(v);
        end
    end
end