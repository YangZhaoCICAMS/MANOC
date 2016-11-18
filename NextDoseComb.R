get.next.manoc<-function(pos.model,j_curr,k_curr,alpha,eta)
{
	ndose.A<-nrow(pos.model)
	ndose.B<-ncol(pos.model)
	
	## Stay. ##
	S.mat<-array(NaN,dim=c(1,3))
	D.mat<-array(NaN,dim=c(2,3))
	E.mat<-array(NaN,dim=c(2,3))
	
	p_S<-pos.model[j_curr,k_curr]
	i_S<-1
	S.mat[i_S,1:2]<-c(j_curr,k_curr)
	S.mat[i_S,3]<-pos.model[j_curr,k_curr]

	## De-escalation. ## 
	p_D<-0
	i_D<-0
	if (j_curr>1)
	{
		i_D<-i_D+1
		D.mat[i_D,1:2]<-c(j_curr-1,k_curr)
		D.mat[i_D,3]<-pos.model[j_curr-1,k_curr]
	}
	
	if (k_curr>1)
	{
		i_D<-i_D+1
		D.mat[i_D,1:2]<-c(j_curr,k_curr-1)
		D.mat[i_D,3]<-pos.model[j_curr,k_curr-1]
	}
	
	if (i_D>0)
	{
		ind_D<-which(D.mat[1:i_D,3]==max(D.mat[1:i_D,3]))
		if (length(ind_D)>1)
		{
			ind_D<-sample(ind_D,size=1,prob=rep(1,length(ind_D)))
		}
		
		D.mat<-D.mat[ind_D,]
		dim(D.mat)<-c(1,3)
		p_D<-D.mat[1,3]
	}
	
	## Escalation. ##
	p_E<-0
	i_E<-0
	if (j_curr<ndose.A)
	{
		i_E<-i_E+1
		E.mat[i_E,1:2]<-c(j_curr+1,k_curr)
		E.mat[i_E,3]<-pos.model[j_curr+1,k_curr]
	}
	
	if (k_curr<ndose.B)
	{
		i_E<-i_E+1
		E.mat[i_E,1:2]<-c(j_curr,k_curr+1)
		E.mat[i_E,3]<-pos.model[j_curr,k_curr+1]
	}
	
	if (i_E>0)
	{
		ind_E<-which(E.mat[1:i_E,3]==max(E.mat[1:i_E,3]))
		if (length(ind_E)>1)
		{
			ind_E<-sample(ind_E,size=1,prob=rep(1,length(ind_E)))
		}
		
		E.mat<-E.mat[ind_E,]
		dim(E.mat)<-c(1,3)
		p_E<-E.mat[1,3]
	}
	
	pp_D<-p_D/(p_D+p_S+p_E)
	pp_S<-p_S/(p_D+p_S+p_E)
	pp_E<-p_E/(p_D+p_S+p_E)
	
	if ((i_D>0)&&(i_E>0))
	{
		p.rgn<-c(pp_D,pp_S,pp_E)
		
		# Dose switching rule. #
		if (max(p.rgn)>=eta)
		{
			action<-which(p.rgn==max(p.rgn))
		}
		else 
		if (max(p.rgn)<eta)
		{
			l<-rep(0,3)
			for (u in 1:3)
			{
				for (v in 1:3)
				{
					l[u]<-l[u]+p.rgn[v]*(alpha*(v-u)*(v>u)+(1-alpha)*(u-v)*(u>v))
				}
			}
			action<-which(l==min(l))
			if (length(action)>1) 
			{
				action<-sample(action,size=1,prob=rep(1,length(action)))
			} 
		}
		
		if (action==1)
		{
            # De-escalation.
			next.dose<-D.mat[1,1:2]
		}
		if (action==2)
		{
			# Retainment.
			next.dose<-S.mat[1,1:2]
		}
		if (action==3)
		{
			# Escalation.
			next.dose<-E.mat[1,1:2]
		}
	}
	
	if (i_D==0)
	{
		p.rgn<-c(pp_S,pp_E)
		
		# Dose switching rule. #
		if (max(p.rgn)>=eta)
		{
			action<-which(p.rgn==max(p.rgn))+1
		}
		else
		if (max(p.rgn)<eta)
		{		
			l<-rep(0,2)
			for (u in 1:2)
			{
				for (v in 1:2)
				{
					l[u]<-l[u]+p.rgn[v]*(alpha*(v-u)*(v>u)+(1-alpha)*(u-v)*(u>v))
				}
			}
			action<-which(l==min(l))+1 # For consistency, action==2-->Retainment. action==3-->Escalation.
		
			if (length(action)>1) 
			{
				action<-sample(action,size=1,prob=rep(1,length(action)))
			} 
		}
		
		if (action==2)
		{
			# Retainment.
			next.dose<-S.mat[1,1:2]
		}
		if (action==3)
		{
			# Escalation.
			next.dose<-E.mat[1,1:2]
		}
	}
	
	if (i_E==0)
	{		
		p.rgn<-c(pp_D,pp_S)
		
		if (max(p.rgn)>=eta)
		{
			action<-which(p.rgn==max(p.rgn))
		}
		else
		if (max(p.rgn)<eta)
		{
			l<-rep(0,2)
			for (u in 1:2)
			{
				for (v in 1:2)
				{
					l[u]<-l[u]+p.rgn[v]*(alpha*(v-u)*(v>u)+(1-alpha)*(u-v)*(u>v))
				}
			}
			action<-which(l==min(l)) # action==1-->De-escalation. action==2-->Stay. action==3-->Escalation.
		
			if (length(action)>1) 
			{
				action<-sample(action,size=1,prob=rep(1,length(action)))
			} 
		}
		if (action==1)
		{
			# De-escalation.
			next.dose<-D.mat[1,1:2]
        }
		if (action==2)
		{
			# Retainment.
			next.dose<-S.mat[1,1:2]
		}
	}
	return(list(next.dose=next.dose))
}
