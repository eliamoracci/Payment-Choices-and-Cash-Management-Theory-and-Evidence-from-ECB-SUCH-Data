function dist=DistanceSMM(param,phi,N,T,truemom)
    beta=param(1);
    eta_h=param(2);
    eta_l=param(3);
    R=param(4);
    sigma=param(5);
    kappa=param(6);
    mu_unc=param(7);
    sigma_unc=param(8);
    simom=GenSimulatedMoments(beta,eta_h,eta_l,R,sigma,kappa,phi,mu_unc,sigma_unc,N,T);
    simom=cell2mat(struct2cell(simom));
    dist=((simom-truemom)./truemom)'*eye(10)*((simom-truemom)./truemom);
end