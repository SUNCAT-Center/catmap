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
        r'(TS|initial_state|final_state|BEP|descriptors)(?:\(([^\(\)]*)\))?(?:\:(?:([\d\w\-]+)|(\[(?:[\d\.]*\,?){1,2}\])))?',
        ['mode','species_list','parameter_key','parameter_list']]

""" Sorry for this regex... I already hate myself for making it! Trying to explain here, but no guarantees
        r'(TS|initial_state|final_state|BEP) <- group1 = mode
        (?: <- terrible python regex notation for open parentheses (fake open parentheses)
        \( <- real open parentheses 
        ([^\(\)]*) <- group2 = species_list
        \) <- real close parentheses
        ) <- end of group
        ? <- I don't know
        (?: <- 
        \: <-colon separates species_list from parameter_list
        (?:([\d\w\-]+)|(\[(?:[\d\.]*\,?){1,2}\])))?' <- this is the parameter list.
"""
