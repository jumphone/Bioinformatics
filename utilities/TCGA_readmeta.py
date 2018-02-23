import sys
import anyjson
fi=open(sys.argv[1])

json=fi.read()

data = anyjson.deserialize(json)

#for one in  data[0]:
#	print one


fo=open(sys.argv[2],'w')
fo.write('#data_type\tfile_name\tfile_id\tsample_type\tsample_id\ttissue_or_organ_of_origin\tprimary_diagnosis\ttumor_stage\tdemographic_id\tproject_id\tcase_id\tstate\tgender\trace\tethnicity\tbirth_year\tdead_year\tdays_to_birth\tdays_to_death\tage_at_diagnosis\tdays_to_recurrence\tdisease_type\n')

def get_key(info,key):
        key="'"+key+"': "
        #key='"'+key+'"'+' '
        info=info[info.find(key)+len(key):]
        if "'," in info:
                info=info[:info.find("',")]
        if "," in info:
                info=info[:info.find(",")]
        if "'}" in info:
                info=info[:info.find("'}")]
        if "}" in info:
                info=info[:info.find("}")]
        if info[:2]=="u'":
                info=info[2:]
        return info.strip().replace(' ','_')




for one in data:
	#if one['file_id']=='cae80695-eeb7-4047-b250-7da7b81ee714':
                info=str(one)
                file_name=one['file_name'] #.replace(' ','_')
                file_id=one['file_id']
                data_type=one['data_type'].replace(' ','_')
                sample_type = get_key(info,"sample_type")
                sample_id = get_key(info,"sample_id")
                tumor_stage=get_key(info,"tumor_stage")
                tissue_or_organ_of_origin=get_key(info,"tissue_or_organ_of_origin")
                primary_diagnosis = get_key(info,"primary_diagnosis")
                demographic_id = get_key(info,"demographic_id")
                p_id=get_key(info,"project_id")
                case_id=get_key(info,"case_id")

                days_to_recurrence=get_key(info,"days_to_recurrence")
                disease_type=get_key(info,"disease_type").replace(' ','_')
                days_to_birth=get_key(info,"days_to_birth")
                gender=get_key(info,"gender")
                state=get_key(info,"state")
                race=get_key(info,"race")
                ethnicity=get_key(info,"ethnicity")
                year_of_death=get_key(info,"year_of_death")
                year_of_birth=get_key(info,"year_of_birth")
                age_at_diagnosis=get_key(info,"age_at_diagnosis")
                #print year_of_death
                #break
                days_to_death=get_key(info,"days_to_death")
                
                fo.write(data_type+'\t'+file_name+'\t'+file_id+'\t'+sample_type+'\t'+sample_id+'\t'+tissue_or_organ_of_origin+'\t'+primary_diagnosis+'\t'+tumor_stage+'\t'+demographic_id+"\t"+p_id+'\t'+case_id+'\t'+state+'\t'+gender+'\t'+race+'\t'+ethnicity+'\t'+year_of_birth+'\t'+year_of_death+'\t'+days_to_birth+'\t'+days_to_death+'\t'+age_at_diagnosis+'\t'+days_to_recurrence+'\t'+disease_type+'\n')	





