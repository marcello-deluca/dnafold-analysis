%   in:
%       COORDS(1xN cell of trajectory frames)
%       FrameRange(Vector of frame numbers to render)
%       StrandNums(vector of particles to which each strand belongs)
%   out:
%       trajectory 

function [] = makeTrajectoryMoviednafold(COORDS,BOXSZ,FrameRange,StrandNums, nstrands, OUTFNAME)
    assert(length(COORDS)~=0,"fed zero-length input coords into movie generator");
    sz = size(COORDS{1});
    writerObj = VideoWriter(OUTFNAME,'MPEG-4');
    writerObj.Quality=100;
    writerObj.FrameRate=24;
    %writerObj.LosslessCompression=true;
    open(writerObj);
    LSIM = length(COORDS);
    if FrameRange==-1
        FrameRange=1:1:LSIM;
    end
    for i = FrameRange
        CURRFRAME = COORDS{i};
        meancrd = mean(CURRFRAME);
        CURR = CURRFRAME - meancrd;
        %CURR = CURRFRAME;
       % color = "black";
       % tubeplot(transpose(CURRFRAME),1,64,[],color);
       %'#2C514C'
       % #057AF0'
       % '#114f11'
       % #3c2ff5
       % '#1f2f8c' FINAL FOR PUBLICATION
       % A00000
       % 4FAD5B
       %CC0000
        hold off
        MakeBeautifiedFrame(CURR(1:512,:),1,'#CC0000');
        %hold off
        %MakeBeautifiedFrame(CURR(1:512,:),1,['#0000ff']);
        hold on
        if StrandNums==-1
            %no other strands
        else
            for q=2:nstrands
                indices=[];
                for j=1:length(StrandNums)
                    if StrandNums(j)==q
                        indices=[indices,j];
                    end
                end
                hold on
                coords = CURR(indices,:);%-meancrd;
                coords = coords-round(coords./BOXSZ(i),0).*BOXSZ(i);
                if (isempty(nonzeros(abs(coords)>BOXSZ(i)*.4)))
                    MakeBeautifiedFrame(coords,1,['#888888']);
                end
            end
        end
        camlight;
        view(-62,5);
        camlight;
        view(-45,5);
        camlight;
        view(-62,5);
        camlight;
        view(88,62);
        camlight;
        ax = gca;
        ax.FontSize = 28;
        
        tstring = sprintf("%.1f ms",i * 10000*10*.01E-9 * 3314 * 1000);
        %RGB = insertText(I,position,tstring);
        %title(tstring);
        txt = text(0,100,65,tstring);
        txt.FontSize=28;
        PLTBOX = 100;
        xlim([-PLTBOX PLTBOX]);
        ylim([-PLTBOX PLTBOX]);
        zlim([-PLTBOX PLTBOX]);
        zoom(1.6);
        %axis equal
        %view(-32,10);
        view(-1,-2)
        fprintf("%i\n",i);
        set(gcf,'color','w');
        frame = getframe(1);
        %im = frame2im(frame);
        %[imind,cm] = rgb2ind(im,256);
        writeVideo(writerObj,frame);
%         if i==FrameRange(1)
%              imwrite(imind,cm,OUTFNAME,'gif', 'Loopcount',inf,'DelayTime',.0416667);
%         else
%              imwrite(imind,cm,OUTFNAME,'gif','WriteMode','append','DelayTime',.0416667);
%         end
    end
    close(writerObj);
end
