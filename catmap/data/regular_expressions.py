regular_expressions = {}

#regular_expressions[key] = [regex,group_name_list]

regular_expressions['site_definition'] = [r'(\d*)(\*)(?:_(\w+))?',
        ['stoichiometry','name','site']]

regular_expressions['species_definition'] = [
        r'(\d*)((?:[^\_\*\+\<\>]+|\*))(?:(?:\*?_(\w+))?|\*)',
        ['stoichiometry','name','site']]

regular_expressions['species_separator'] = ['(?:\A|\s*\+\s*|\s+)',[]]

regular_expressions['initial_transition_final_states'] = [
        r'([^\<\>]*)(?:\<?\-\>)(?:([^\<\>]*)(?:\<?\-\>))?([^\<\>]*)',
        ['initial_state','transition_state','final_state']]

regular_expressions['transition_state_scaling_constraint'] = [
        r'(TS|initial_state|final_state|BEP)(?:\(([^\(\)]*)\))?(?:\:(?:([\d\w\-]+)|(\[(?:[\d\.]*\,?){1,2}\])))?',
        ['mode','species_list','parameter_key','parameter_list']]
