% Figure 2: Triangular phases at 8different 
clear;

CALCON = 0;
filename1 = '../data/old/figureS2ATriMesh';
filename2 = '../data/old/figureS2BTriMesh';
filename3 = '../data/old/figureS2CTriMesh';
filename4 = '../data/old/figureS2DTriMesh';
filenames = {filename1, filename2, filename3, filename4};

if CALCON
    % Figure 2A-D
    cnt = 1;
    for N = [1e0, 1e3]
        for CHI = [0, 0.8]
            filename = filenames{cnt};
            
            [FAVV, PHIPV, EIGV, EIG, KSV] = calcphase_wsolvent(N, CHI, 1, 0);
            save(filename, 'FAVV', 'PHIPV', 'EIGV', 'EIG', 'KSV', 'N');
            cnt = cnt+1;
        end
    end
    
else
    for ii = 1:4
        filename = filenames{ii};
        load(filename);
        
        figure;plotphase_wsolvent(FAVV, PHIPV, EIGV, EIG, KSV, N);
    end
end