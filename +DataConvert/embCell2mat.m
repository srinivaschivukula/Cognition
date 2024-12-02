function out=embCell2mat(in)

if ndims(in{1})==2
    
    if any(size(in{1})==1)
        for i=1:length(in{1})
            out(:,i)=squeeze(cell2mat(cellfun(@(x)x(i), in,'UniformOutput',false)));
        end
        
    else
        for i=1:length(in)
            out(:,:,i)=in{i};
        end
    end
end