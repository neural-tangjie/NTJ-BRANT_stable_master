function Graph_Connect(x,y)
%����: xΪһ��int����,����n=5,x=[1,2,3,4,5]
%      yΪ��άint����,����y=[2,3,1]����y=[5,2,1]
%  ���:��ͼ������
%��1��ÿ��x�е�������һ��СԲȦȦ����,Ȼ����ȷֲ���һ�������ԲȦ�ϣ�ԲȦ����Ҫ���֣� 
%��2������y�е�ֵ��һ������ͷ���߽��������ڵ�Ȧ������������y��ǰ��������Ӧx�еĵ㣬����������Ϊ1������y��ǰ��������˳����Ӵ���ͷ�����ߣ�����������Ϊ�㣬���Ѿ��е�����ȥ����
%          ��y=[2,3,1]����2-----��3,  ���һ����2��3�ļ�ͷ����
%          ��y=[5,2,1]����5-----��2   ���һ����5��2�ļ�ͷ����
%          ��y=[2,5,1]����5<-----��2 ���һ����2��5�ļ�ͷ���� 
%          ��y=[5,2,0]����5<----- 2   ɾ��һ����5��2�ļ�ͷ����

%���x�Ƿ�Ϊ��һ��������
figure(1);hold on;
if (size(x,1)>size(x,2))
    x=x';
end
n=size(x,2);%n:��ʾx�����Ȼ�������ĳ���

%��x�е�ÿ���������ȵķֲ���һ�������ԲȦ��
vir_r=10;%����ԲȦ�İ뾶
r=0.5;%Ȧ�����ֵ�СԲȦ�İ뾶

dita=linspace(0,2*pi,n+1);

sita=0:pi/20:2*pi;
% figure(1);
extend=1;% ����ͼ��������һ��
axis([-vir_r-extend,vir_r+extend,-vir_r-extend,vir_r+extend]);
hold on;

for i=1:n 
    str=int2str(x(i));
    %vir_r*cos(dita(i))-------------СԲȦԲ�ĵĺ�����
    %vir_r*sin(dita(i)--------------СԲ--------������
    text(vir_r*cos(dita(i)),vir_r*sin(dita(i)),str);%��СԲȦ��д����
    %��СԲȦ��Բ�Ĵ�����һ���뾶Ϊr��Բ
    plot(vir_r*cos(dita(i))+r*cos(sita),vir_r*sin(dita(i))+r*sin(sita)); 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����yֵ���Ƽ�ͷ
%���y�Ƿ�Ϊһ��������,�����򷵻�
if(size(y,1)>size(y,2))
    return;
end

lable=y(3);

%���߳�����Ҫ˼·: ���߳������yֵ,���Ӵ���㵽�յ��СԲԲ��,���һ��Ƽ�ͷ
%
%
arrowlength=r;%��ͷ���ڵȱ������εĸߵĳ���
if(lable==1)%���߳���
    startpoint_x=vir_r*cos(dita(y(1)));
    startpoint_y=vir_r*sin(dita(y(1))); 
    
    %��Ϊ��ͷ��Ҫ�˵�ԲȦ����,��ʱΪ��ʱ�յ�
    tempendpoint_x=vir_r*cos(dita(y(2)));
    tempendpoint_y=vir_r*sin(dita(y(2)));
    
    
    k=(tempendpoint_y-startpoint_y)/(tempendpoint_x-startpoint_x);%б��
    
    if(startpoint_x<tempendpoint_x)%�������һ���
        
        if(k==0)
            endpoint_x=tempendpoint_x-r;
            endpoint_y=tempendpoint_y;
            plot([startpoint_x endpoint_x],[startpoint_y endpoint_y]); 
            
            plot([endpoint_x endpoint_x-arrowlength*0.866],[endpoint_y endpoint_y+arrowlength*0.5]);
            plot([endpoint_x endpoint_x-arrowlength*0.866],[endpoint_y endpoint_y-arrowlength*0.5]);
        else
             angle=atan(k);
            %����ˮƽ���ʱ,�յ���˺����
            endpoint_x=tempendpoint_x-r*cos(angle);
            endpoint_y=startpoint_y+k*(endpoint_x-startpoint_x);
            
            plot([startpoint_x endpoint_x],[startpoint_y endpoint_y]); 
            
            %���ݽ������ι�ϵ,���Ƽ�ͷ���ز���
            uparrow_x=endpoint_x-arrowlength*cos(angle)-arrowlength*1.732/3*sin(angle);
            uparrow_y=startpoint_y+k*(endpoint_x-arrowlength*cos(angle)-startpoint_x)+1/k*1.732/3*arrowlength*sin(angle);            
            plot([endpoint_x uparrow_x],[endpoint_y uparrow_y]);
            %plot([endpoint_x-r endpoint_x-r-arrowlength*0.866],[endpoint_y endpoint_y+arrowlength*0.5]);
            
            %���ݽ������ι�ϵ,���Ƽ�ͷ���ز���
            downarrow_x=endpoint_x-arrowlength*cos(angle)+arrowlength*1.732/3*sin(angle);
            downarrow_y=startpoint_y+k*(endpoint_x-arrowlength*cos(angle)-startpoint_x)-1/k*1.732/3*arrowlength*sin(angle);
            plot([endpoint_x downarrow_x],[endpoint_y downarrow_y],'r');          
        end
    end
    
    if(startpoint_x>tempendpoint_x)%�����������
        if(k==0)
            endpoint_x=tempendpoint_x+r;
            endpoint_y=tempendpoint_y;
            plot([startpoint_x endpoint_x],[startpoint_y endpoint_y]); 
            
            plot([endpoint_x endpoint_x+arrowlength*0.866],[endpoint_y endpoint_y+arrowlength*0.5]);
            plot([endpoint_x endpoint_x+arrowlength*0.866],[endpoint_y endpoint_y-arrowlength*0.5]);
        else
             angle=atan(k)+pi;%�����к�������ֵ��[-0.5*pi,0.5*pi]֮��,+pi�ָ���ԭ���ĽǶ�
            %����ˮƽ���ʱ,�յ��ǰ�����
            endpoint_x=tempendpoint_x-r*cos(angle);
            endpoint_y=startpoint_y+k*(endpoint_x-startpoint_x);
            
            plot([startpoint_x endpoint_x],[startpoint_y endpoint_y]); 
            
            %���ݽ������ι�ϵ,���Ƽ�ͷ���ز���
            uparrow_x=endpoint_x-arrowlength*cos(angle)-arrowlength*1.732/3*sin(angle);
            uparrow_y=startpoint_y+k*(endpoint_x-arrowlength*cos(angle)-startpoint_x)+1/k*1.732/3*arrowlength*sin(angle);            
            plot([endpoint_x uparrow_x],[endpoint_y uparrow_y]);
            %plot([endpoint_x-r endpoint_x-r-arrowlength*0.866],[endpoint_y endpoint_y+arrowlength*0.5]);
            
            %���ݽ������ι�ϵ,���Ƽ�ͷ���ز���
            downarrow_x=endpoint_x-arrowlength*cos(angle)+arrowlength*1.732/3*sin(angle);
            downarrow_y=startpoint_y+k*(endpoint_x-arrowlength*cos(angle)-startpoint_x)-1/k*1.732/3*arrowlength*sin(angle);
            plot([endpoint_x downarrow_x],[endpoint_y downarrow_y],'r');         
        end
    end
end

if(lable==0)%����������ͷ�ĳ���,��ȫ����������ͼ����,ֻ����һ�λ��Ƶ���ͬ������ɫһ�µĴ���ͷֱ��
    startpoint_x=vir_r*cos(dita(y(1)));
    startpoint_y=vir_r*sin(dita(y(1))); 
    
    %��Ϊ��ͷ��Ҫ�˵�ԲȦ����,��ʱΪ��ʱ�յ�
    tempendpoint_x=vir_r*cos(dita(y(2)));
    tempendpoint_y=vir_r*sin(dita(y(2)));
    
    
    k=(tempendpoint_y-startpoint_y)/(tempendpoint_x-startpoint_x);%б��
    
    if(startpoint_x<tempendpoint_x)%�������һ���
        
        if(k==0)
            endpoint_x=tempendpoint_x-r;
            endpoint_y=tempendpoint_y;
            plot([startpoint_x endpoint_x],[startpoint_y endpoint_y],'w'); 
            
            plot([endpoint_x endpoint_x-arrowlength*0.866],[endpoint_y endpoint_y+arrowlength*0.5],'w');
            plot([endpoint_x endpoint_x-arrowlength*0.866],[endpoint_y endpoint_y-arrowlength*0.5],'w');
        else
             angle=atan(k);
            %����ˮƽ���ʱ,�յ���˺����
            endpoint_x=tempendpoint_x-r*cos(angle);
            endpoint_y=startpoint_y+k*(endpoint_x-startpoint_x);
            
            plot([startpoint_x endpoint_x],[startpoint_y endpoint_y],'w'); 
            
            %���ݽ������ι�ϵ,���Ƽ�ͷ���ز���
            uparrow_x=endpoint_x-arrowlength*cos(angle)-arrowlength*1.732/3*sin(angle);
            uparrow_y=startpoint_y+k*(endpoint_x-arrowlength*cos(angle)-startpoint_x)+1/k*1.732/3*arrowlength*sin(angle);            
            plot([endpoint_x uparrow_x],[endpoint_y uparrow_y],'w');
            %plot([endpoint_x-r endpoint_x-r-arrowlength*0.866],[endpoint_y endpoint_y+arrowlength*0.5]);
            
            %���ݽ������ι�ϵ,���Ƽ�ͷ���ز���
            downarrow_x=endpoint_x-arrowlength*cos(angle)+arrowlength*1.732/3*sin(angle);
            downarrow_y=startpoint_y+k*(endpoint_x-arrowlength*cos(angle)-startpoint_x)-1/k*1.732/3*arrowlength*sin(angle);
            plot([endpoint_x downarrow_x],[endpoint_y downarrow_y],'w');          
        end
    end
    
    if(startpoint_x>tempendpoint_x)%�����������
        if(k==0)
            endpoint_x=tempendpoint_x+r;
            endpoint_y=tempendpoint_y;
            plot([startpoint_x endpoint_x],[startpoint_y endpoint_y],'w'); 
            
            plot([endpoint_x endpoint_x+arrowlength*0.866],[endpoint_y endpoint_y+arrowlength*0.5],'w');
            plot([endpoint_x endpoint_x+arrowlength*0.866],[endpoint_y endpoint_y-arrowlength*0.5],'w');
        else
             angle=atan(k)+pi;%�����к�������ֵ��[-0.5*pi,0.5*pi]֮��,+pi�ָ���ԭ���ĽǶ�
            %����ˮƽ���ʱ,�յ��ǰ�����
            endpoint_x=tempendpoint_x-r*cos(angle);
            endpoint_y=startpoint_y+k*(endpoint_x-startpoint_x);
            
            plot([startpoint_x endpoint_x],[startpoint_y endpoint_y],'w'); 
            
            %���ݽ������ι�ϵ,���Ƽ�ͷ���ز���
            uparrow_x=endpoint_x-arrowlength*cos(angle)-arrowlength*1.732/3*sin(angle);
            uparrow_y=startpoint_y+k*(endpoint_x-arrowlength*cos(angle)-startpoint_x)+1/k*1.732/3*arrowlength*sin(angle);            
            plot([endpoint_x uparrow_x],[endpoint_y uparrow_y],'w');
            %plot([endpoint_x-r endpoint_x-r-arrowlength*0.866],[endpoint_y endpoint_y+arrowlength*0.5]);
            
            %���ݽ������ι�ϵ,���Ƽ�ͷ���ز���
            downarrow_x=endpoint_x-arrowlength*cos(angle)+arrowlength*1.732/3*sin(angle);
            downarrow_y=startpoint_y+k*(endpoint_x-arrowlength*cos(angle)-startpoint_x)-1/k*1.732/3*arrowlength*sin(angle);
            plot([endpoint_x downarrow_x],[endpoint_y downarrow_y],'w');         
        end
    end
end
hold off;