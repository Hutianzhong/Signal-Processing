function y = operator(f,str,length,high)

refSE = strel('line',length,0);
    switch str
        case '+'
            for le = 1 : high
                for lee = 1 : le
                f(le,:) = imdilate(f(le,:),refSE);
                end
            end
            
        case '-'
            for le = 1 : high
                for lee = 1 : le
                f(le,:) = imerode(f(le,:),refSE);
                end
            end
    end
    
    y = f;
end