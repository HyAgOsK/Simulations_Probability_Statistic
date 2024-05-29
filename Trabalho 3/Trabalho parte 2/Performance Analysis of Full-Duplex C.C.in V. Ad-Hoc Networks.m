
% Simula��o desenvolvida para a disciplina TP547 - Princ�pios de Simula��o de
% Sistemas de Comunica��o do curso de Mestrado do Inatel
% Igor Gon�alves de Souza 931
% Hyago Vieira Lemes Barbosa Silva 922
% Performance Analysis of Full-Duplex Cooperative Communication in Vehicular
% Ad-Hoc Networks
clc, clear all, close all, fprintf('\n')
warning('off', 'all');
pkg load signal

% Vari�ncia do ru�do e bits por uso do canal.
% Coeficiente de Perda de Percurso e Pot�ncia m�dia.
N0 = 1; R = 3;	 alpha = 4; y_rlrl = 10^(-4);
N = 1e6; 		% Amostras Monte Carlo.

% Dist�ncias entre a Fonte (s), o Rel� (rl) e o Destino (d).
d_srl = 0.5; d_rld = 0.5; d_sd = 1;

% O canal hij, com i {s, rl} e j {rl, d} � hij = sqrt(yij) * desvanecimento
% com distribui��o Nakagami-m modelado pelo par�metro de desvanecimento
% mij e pot�ncia m�dia yij.
% A simula��o considera o Caso 1: m_sd = 0.5 e m_srl = m_rlrl = m_rld = 1.
h_srl = (max(0, sqrt(d_srl)^(-alpha).*randraw('nakagami', 1, [1 N]))).^2;
h_rld = (max(0, sqrt(d_rld)^(-alpha).*randraw('nakagami', 1, [1 N]))).^2;
h_sd = (max(0, sqrt(d_sd)^(-alpha).*randraw('nakagami', 0.5, [1 N]))).^2;
h_rlrl = (max(0, sqrt(y_rlrl).*randraw('nakagami', 1, [1 N]))).^2;

O_VJD = {}; O_VDH = {}; O_VHD = {}; P = -5:20;
for p = -5:20
		Ps = db2pow(p); P_rl = Ps;

		% O esquema VJD conta com a ajuda de rel� FD de forma que o link direto
		% seja visto como informa��o �til em vez de interfer�ncia no destino.
		Isrl_VJD = log2(1 + (Ps.*h_srl)./(P_rl.*h_rlrl + N0)); 		% Equa��o (5).
		Isd_VJD = log2(1 + (Ps.*h_sd)./N0);														% Equa��o (6).
		Irld_VJD = log2(1 + (Ps.*h_sd + P_rl.*h_rld)./N0);				% Equa��o (7).
		O_VJD{end + 1} = sum(max(Isd_VJD, min(Isrl_VJD, Irld_VJD)) < R)./N;

		% No esquema VDH o link direto sd � visto como interfer�ncia no destino.
		% A informa��o m�tua do link sr � escrita como em (5), ou seja,
		% Isrl_VDH = Isrl_VJD.
		Isrl_VDH = Isrl_VJD;
		Irld_VDH = log2(1 + (P_rl.*h_rld)./(Ps.*h_sd + N0));	% Equa��o (10).
		O_VDH{end + 1} = sum(min(Isrl_VJD,  Irld_VDH) < R)./N;

		% No esquema VHD, a transmiss�o ocorre em dois intervalos de tempo.
		% No primeiro, a fonte transmite sua mensagem para o rel� selecionado e
		% para o destino, enquanto no segundo intervalo de tempo o rel� retransmite
		% a mensagem de origem se for corretamente decodificada e solicitada pelo
		% destino.
		Isd_VHD = 0.5*log2(1 + (Ps.*h_sd)./N0);												% Equa��o (13).
		Isrl_VHD = 0.5*log2(1 + (Ps.*h_srl)./N0);												% Equa��o (13).
		Irld_VHD = 0.5*log2(1 + (Ps.*h_sd + P_rl.*h_rld)./N0);	% Equa��o (14).
		O_VHD{end + 1} = sum(max(Isd_VHD, min(Isrl_VHD, Irld_VHD)) < R)./N;
end

% O Gr�fico da Figura 3 apresenta a probabilidade de falha em fun��o da
% pot�ncia de transmiss�o P, em que d_srl = d_rld = 0.5, d_sd = 1,
% y_rlrl =10^(-4) e R = 3.
line = {'linewidth', 1.0}; style = {'k-.', 'r-', 'g--'};
semilogy(P, cell2mat(O_VDH), style{1}, line{:},  P, cell2mat(O_VHD), style{2},  line{:},
						   P, cell2mat(O_VJD), style{3},  line{:}), grid on
title({'Performance Analysis of Full-Duplex Cooperative Communication',  'in Vehicular Ad-Hoc Networks'})
xlabel('\fontsize{14}{P [db]}'), ylabel('\fontsize{14}{\it{Outage Probability}}')
legend('VDH', 'VHD', 'VJD', 'Location', 'southwest')
print -dpng AdHoc_Networks.png

