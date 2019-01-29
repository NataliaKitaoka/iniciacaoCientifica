//PCA das imagens de colônias de S. cerevisiae irradiadas com 2Gy e do controle. 
clear, clear all, clearglobal()
//stacksize(2e7)

//------------- compressao por Wavelet -------------------

n= 34;// numero de dados totais
n1 = 12; //numero de amostras da classe 1 = controle 
n2 = 8; //numero de amostras da classe 2.1 = irradiadas placa 1 
n3 = 14; //numero de amostras da classe 2.2 = irradiadas placa 2 
n4 = n2+n3; //total de irradiadas = 22
x1 = dir('C:\Users\natal\Desktop\IC\experimentosImagens\15.01.19 - irradiada\controle'); //lista dos arquivos na pasta de imagens
x2 = dir('C:\Users\natal\Desktop\IC\experimentosImagens\15.01.19 - irradiada\placa1_2Gy');
x3 = dir('C:\Users\natal\Desktop\IC\experimentosImagens\15.01.19 - irradiada\placa2_2Gy');

//leitura e compressao de todos os arquivos da pasta controle
for i = 1:n1
    c = sci2exp(i)// converte constante em string
    im0 = imread('C:\Users\natal\Desktop\IC\experimentosImagens\15.01.19 - irradiada\controle\' + x1(i)); //para cada teste criar uma pasta com todas as amostras que se queira
    im1 = rgb2gray(im0);
    im2=im2double(im1);
    [C,S]=wavedec2(im2,3,'bior1.1');
    cA=appcoef2(C,S,'bior1.1',3);
    a=wkeep(cA,[size(cA)]); //mesma coisa que cA
    a=(a-min(a))/max(a-min(a)); 
    m1 =sum(a, 'r');
    m1=(m1-min(m1))/max(m1-min(m1));  
    save('C:\Users\natal\Desktop\IC\experimentosImagens\15.01.19 - irradiada\dados.sce'+c+'.dat');//salvar em outra pasta criada (pasta de dados)
    
end

//leitura e compressao de todos os arquivos da pasta irradiadas placa 1
for j = 1:n2
    c = sci2exp(j+i)// converte constante em string
    im0 = imread('C:\Users\natal\Desktop\IC\experimentosImagens\15.01.19 - irradiada\placa1_2Gy\' + x2(j)); //para cada teste criar uma pasta com todas as amostras que se queira
    im1 = rgb2gray(im0);
    im2=im2double(im1);
    [C,S]=wavedec2(im2,3,'bior1.1');
    cA=appcoef2(C,S,'bior1.1',3);
    a=wkeep(cA,[size(cA)]); //mesma coisa que cA
    a=(a-min(a))/max(a-min(a)); 
    m1 =sum(a, 'r');
    m1=(m1-min(m1))/max(m1-min(m1));  
    save('C:\Users\natal\Desktop\IC\experimentosImagens\15.01.19 - irradiada\dados.sce'+c+'.dat');//salvar em outra pasta criada (pasta de dados)
    
end

//leitura e compressao de todos os arquivos da pasta irradiadas placa 2
for l = 1:n3
    c = sci2exp(l+j)// converte constante em string
    im0 = imread('C:\Users\natal\Desktop\IC\experimentosImagens\15.01.19 - irradiada\placa2_2Gy\' + x3(l)); //para cada teste criar uma pasta com todas as amostras que se queira
    im1 = rgb2gray(im0);
    im2=im2double(im1);
    [C,S]=wavedec2(im2,3,'bior1.1');
    cA=appcoef2(C,S,'bior1.1',3);
    a=wkeep(cA,[size(cA)]); //mesma coisa que cA
    a=(a-min(a))/max(a-min(a)); 
    m1 =sum(a, 'r');
    m1=(m1-min(m1))/max(m1-min(m1));  
    save('C:\Users\natal\Desktop\IC\experimentosImagens\15.01.19 - irradiada\dados.sce'+c+'.dat');//salvar em outra pasta criada (pasta de dados)
    
end
// --------------------- PCA -------------------------

X = [];

for i = 1:n
    c = sci2exp(i)// converte constante em string
    load('C:\Users\natal\Desktop\49IMAGENScomprimidas\'+c+'.dat'); //usar mesmo diretório da pasta de dados
    X = [X m1'] //X = matriz de dados
end

m = size(X);
m = m(1);

//Calculando a média de cada linha de X
for i=1:m
    Xm(i)=mean(X(i,:)); // retorna um vetor coluna
    Xstd(i)=stdev(X(i,:));
end


//subtrai a media de cada coluna
for j=1:n//numero de amostras
    Xmm(:,j)=X(:,j)-Xm;
end
  
Xmm;    
    
for i=1:m
    Xmmm(i,:)=Xmm(i,:)/Xstd(i,:); // divide cada coluna pelo desvio padrao
end

Xmmm;

//PCAScatterPlot(Xmmm', 'r', '.') CALCULA A PCA

//Calculando a matrix de covariância

Xc=(1/n)*(Xmmm*Xmmm'); // dividir pelo número de amostras


//Encontrando os autovetores (componentes principais) e autovalores

[PC1,Av1]=spec(Xc);

//Extrair a diagonal da matriz dos autovalores como um vetor

Av1=diag(Av1);

//Organizar as variancias em ordem decrescente
[J1,I1]=gsort(Av1);

//Organizar os autovetores em ordem decrescente de representatividade dos dados
PC1=PC1(:,I1);

//Projetar os dados originais nas componentes principais
Sig1=PC1'*Xmmm;
//Sig1 = abs(Sig1)
//Sig1 = real(Sig1)

//Plotar os dados nas componentes PC1 e PC1
figure(1),plot(Sig1(1,:),Sig1(1,:),'r^');

//Plotar os dados nas componentes PC1 e PC2
figure(2),plot(Sig1(1,1:n1),Sig1(2,1:n1),'r^');
figure(2),plot(Sig1(1,n1+1:n4),Sig1(2,n1+1:n4),'b^');
//figure(2),plot(Sig1(1,26:39),Sig1(2,26:39),'g^');
//figure(2),plot(Sig1(1,40:49),Sig1(2,40:49),'y^');

//Plotar os dados nas componentes PC1 e PC3
figure(3),plot(Sig1(1,:),Sig1(3,:),'r^');

//Plotar os dados nas PC 1,2 e 3
////B. subtilis (1- 13); S. aureus (14-25 imagens); S.cerevisiae (26-39); E. coli (40-49).
//r= B. subtilis, b=S. aureus, g=S. cerevisiae e y=E.coli
//figure(4),plot3d(Sig1(1,1:13),Sig1(2,1:13),Sig1(3,1:13),'r^');
//figure(4),plot3d(Sig1(1,14:25),Sig1(2,14:25),Sig1(3,14:13),'b^');
//figure(4),plot3d(Sig1(1,26:39),Sig1(2,26:39),Sig1(3,26:39),'g^');
//figure(4),plot3d(Sig1(1,40:49),Sig1(2,40:49),Sig1(3,40:49),'y^');

figure(5), plot(J1/sum(J1));
