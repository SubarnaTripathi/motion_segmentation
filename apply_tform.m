function Iw=apply_tform(I,T,tformtype)
    % function Iw=apply_affine(I,A,t)
    %  apply the A,t affine transformation to I

    switch tformtype
    case 'affine',
      Tf = maketform('affine',T');
      Iw=imtransform(I,Tf,'xdata',[1 size(I,2)],'ydata',[1 size(I,1)], ...
             'fillvalues',NaN);

    case 'projective',
      Tf = maketform('projective',T');
      Iw=imtransform(I,Tf,'xdata',[1 size(I,2)],'ydata',[1 size(I,1)],'fillvalues',NaN);
    otherwise
      Iw=I;
    end
end


