function standard_comp_name=get_comp_name(comp,file_name)
% Gets the sac component name and filename and try to identify e, n and z components
    %% Bidouille pour les noms de fichier spectre, dans le cas des donnees de Marcelo, les noms de fichiers
    %% et les noms de composantes dans les headers sac ne sont pas homogenes.
    %% Hypothese fichiers .2.sac.vel et .5.sac.vel sont N, et fichiers .3.sac.vel et .6.sac.vel sont E
    %% De plus si comp=radial je modifie pour nord et transversal pour est
    
    if strcmp(comp(1:2),'HH') ~= 1
        if strcmp(comp,'vertical') == 1 % may add v, z, Z in the test
            standard_comp_name='z';
        elseif strcmp(comp,'north') == 1 % may add n, N in the test
            standard_comp_name='n';
        elseif strcmp(comp,'east') == 1 % may add e, E in the test
            standard_comp_name='e';
        else
            % the following works with some files send by Marcelo but needs a warning anyway
            % see disp below
            % If component name not recognised, look at filename
            % Assumption the file is blabla.comp_info.sac
            % then gets the comp_info field from filename and compare with "known" codes
            % i.e. 1, 4, z, Z for vertical
            %      2, 5, n, N for north-south
            %      3, 6, e, E for east-west
            pos=strfind(file_name,'.');
            if length(pos) > 1
                test=file_name(pos(length(pos)-1)-1:pos(length(pos)-1)-1);
                if (strcmp(test,'1') == 1 || strcmp(test,'4') == 1 || strcmp(test,'z') == 1 || strcmp(test,'Z') == 1)
                    standard_comp_name='z';
                elseif (strcmp(test,'2') == 1 || strcmp(test,'5') == 1 || strcmp(test,'n') == 1 || strcmp(test,'N') == 1 || strcmp(test,'r') == 1)
                    standard_comp_name='n';
                elseif (strcmp(test,'3') == 1 || strcmp(test,'6') == 1 || strcmp(test,'e') == 1 || strcmp(test,'E') == 1 || strcmp(test,'t') == 1)
                    standard_comp_name='e';
                else
                    standard_comp_name=-99;
                end
                if standard_comp_name ~= -99
                    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                    disp('warning component recognition based on filename')
                    disp(file_name)
                    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                end
            else
                standard_comp_name=-99;
            end
        end
    else
        if strcmp(comp,'HHZ') == 1
            standard_comp_name='z';
        elseif strcmp(comp,'HHN') == 1
            standard_comp_name='n';
        elseif strcmp(comp,'HHE') == 1
            standard_comp_name='e';
        end
    end
endfunction