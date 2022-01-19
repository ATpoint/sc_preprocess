#! /usr/bin/env nextflow 

nextflow.enable.dsl=2

/*********************************************  

    PARAMS VALIDATION VIA SCHEMA.NF

**********************************************/

// Terminal colors:
ANSI_RESET   = "\u001B[0m"
ANSI_RED     = "\u001B[31m"
ANSI_GREEN   = "\u001B[32m"
ANSI_YELLOW  = "\u001B[33m"
DASHEDDOUBLE = "=".multiply(70)
    
// Function for printing consistent error messages:
def ErrorMessenger(base_message='', additional_message=''){
    println("$ANSI_RED" + "$DASHEDDOUBLE")
    println "[VALIDATION ERROR]"
    println base_message
    if(additional_message!='') { println("$additional_message") }
    println("$DASHEDDOUBLE" + "$ANSI_RESET")
}

// Main validation function:
def ValidateParams(){
    
    // Location of the schema file, expected in the main dir:
    def schemafile = new File("$projectDir/schema.nf")

    // [VALIDATION] Check that schema.nf exists
    if(!schemafile.exists()){
        ErrorMessenger("The expected \$baseDir/schema.nf file does not exist!",
                       "=> File 'schema.nf' must be located in the same directory as the main.nf!")
        System.exit(1)
    }

    // Import schema map from schema.nf
    def schema = evaluate(schemafile)

    // Parse keys from params and schema:
    def schema_keys = schema.keySet()
    def params_keys = params.keySet()
    def schema_error = 0

    // [VALIDATION] Each param from "standard" Nextflow must be defined in schema.nf:
    def diff_keys   = params_keys - schema_keys
    if(diff_keys.size()>0){
        ErrorMessenger("These params are set but not defined in schema.nf:",
                       diff_keys.each { k -> "--${k}"})
        schema_error += 1
    }
    
    // Go through each schema entry and validate:
    schema.each { schema_name, entry -> 
    
        // If there is a param with the same name use the 'value' of it instead of the schema value:
        if(params.keySet().contains(schema_name)){
            schema[schema_name]['value'] = params[schema_name]
        }
        
        def schema_value     = entry['value']
        def schema_type      = entry['type']
        def schema_allowed   = entry['allowed']
        def schema_mandatory = entry['mandatory']
        def schema_pattern   = entry['pattern']

        // [VALIDATION] The two keys 'value' and 'type' must be set:
        def keys_diff = ["value", "type"] - entry*.key
        if(keys_diff.size() > 0){
            ErrorMessenger("schema.${schema_name} does not contain the four keys value/type/mandatory/allowed",
                           "=> The schema map must look like: schema.foo = [value:, type:, mandatory:, allowed:]")
            schema_error += 1
            return
        }   

        // [VALIDATION] Only the five allowed keys [value, type, allowed, mandatory, pattern] can be set:
        def allowed_keys = ["value", "type", "mandatory", "allowed", "pattern"]
        def keys_diff2 = entry*.key - allowed_keys
        if(keys_diff2.size() > 0){
            ErrorMessenger("schema.${schema_name} contains non-allowed keys!",
                           "=> Allowed keys are: ${allowed_keys}")
            schema_error += 1
            return
        } 

        // [VALIDATION] If 'mandatory' is true then 'value' must not be empty
        def mandatory_not_empty_error = "schema.${schema_name} is mandatory but not set or empty"
        if(schema_mandatory){
            if(schema_value.toString().trim()=='' || schema_value.toString().trim()==null) {
                ErrorMessenger(mandatory_not_empty_error) 
                schema_error += 1
                return
            } 
        }

        // [VALIDATION] The 'type' must be one of the allowed choices for this key:
        def type_allowed_choices = ['integer', 'float', 'numeric', 'string', 'logical']
        def type_allowed_choices_error = "The 'type' key in schema.${schema_name} must be one of:"
        def type_allowed_choices_valid = type_allowed_choices.contains(schema_type)
        if(!type_allowed_choices_valid){
            ErrorMessenger(type_allowed_choices_error, type_allowed_choices)
            schema_error += 1
            return
        }

        // [VALIDATION] The 'value' must have the correct 'type':
        def value_type_match_error = "schema.${schema_name} is not of type $schema_type"
        if(schema_type=="integer"){
            if((schema_value !instanceof Integer) && (schema_value !instanceof Long)){
                ErrorMessenger(value_type_match_error, "=> You provided: $schema_value")
                schema_error += 1
                return
            }
        }                 
            
        if(schema_type=="float"){
            if((schema_value !instanceof Double) && (schema_value !instanceof Float) && (schema_value !instanceof BigDecimal)){
                ErrorMessenger(value_type_match_error, "=> You provided: $schema_value")
                schema_error += 1
                return
            }                 
        }

        if(schema_type=="numeric"){
            if((schema_value !instanceof Integer) && (schema_value !instanceof Long) &&
               (schema_value !instanceof Double) && (schema_value !instanceof Float) && 
               (schema_value !instanceof BigDecimal)){
                ErrorMessenger(value_type_match_error, "=> You provided: $schema_value")
                schema_error += 1
                return
            }                 
        }                                  

        if(schema_type=="string"){
            if((schema_value != '') && (schema_value !instanceof String) && (schema_value !instanceof GString)){
                ErrorMessenger(value_type_match_error, "=> You provided: $schema_value")
                schema_error += 1
                return
            }                 
        }

        if(schema_type=="logical"){
            if(schema_value !instanceof Boolean){
                ErrorMessenger(value_type_match_error, "=> You provided: $schema_value")
                schema_error += 1
                return
            }                 
        }

        // [VALIDATION] The 'value' must be part of the 'allowed' choices:
        if(schema_allowed!=null){
            allowed_not_contain_value_error = "The 'value' of ${schema_name} is not allowed!"
            if(schema_allowed!=''){
                if(!schema_allowed.contains(schema_value)){
                    ErrorMessenger(allowed_not_contain_value_error,
                                   "=> You provided: ${schema_value}\nAllowed options are:\n${schema_allowed}")
                    schema_error += 1
                    return
                }
            }
        }
        
        // [VALIDATION] schema_allowed choices must be of same 'type' as 'value'
        if(schema_allowed!=null){
            def allowed_type_match_error = "schema.${schema_name} must be one of: \n${schema_allowed}"
            if(schema_allowed!=''){
                def value_class = schema_value.getClass()
                schema_allowed.each { 
                    def current_class = it.getClass()
                    if(current_class != value_class) {
                        ErrorMessenger(allowed_type_match_error) 
                        schema_error += 1
                        return
                    }
                }  
            }
        }

        // [VALIDATION] The 'value' obeys 'pattern':
        if(schema_pattern!=null & schema_type=='string' & schema_value!=''){
            def value_not_obey_pattern_error = "The 'value' of schema.${schema_name} does not match the 'pattern'"
            def value_not_obey_pattern = schema_value ==~ schema_pattern
            if(!value_not_obey_pattern){
                ErrorMessenger(value_not_obey_pattern_error,
                              "=> The expected pattern is ${schema_pattern}")
                schema_error += 1
                return
            }
        }
        
        // All validations successful, now feed the 'value' into the global params to be used in the workflows:
        params[schema_name] = schema_value

    }

    // If there were any validation fails print a summary error message and exit:
    if(schema_error > 0){

        def was_were = schema_error==1 ? "was" : "were"
        def spacer = was_were=="was" ? "      " : "    "
        def xerrors = schema_error==1 ? "error" : "errors"
        println("$ANSI_RED" + "$DASHEDDOUBLE")
        println "||                                                                  ||"
        println "||                                                                  ||"
        println("||          [EXIT ON ERROR] Parameter validation failed!            ||")
        println "||                                                                  ||"
        println("||      There $was_were a total of $schema_error validation $xerrors for schema.nf!${spacer}||")
        println "||                                                                  ||"
        println "||                                                                  ||"
        println("$DASHEDDOUBLE" + "$ANSI_RESET")

        System.exit(1)     

    } else {

        println("${ANSI_GREEN}${DASHEDDOUBLE}")
        println("[INFO] Schema validation successful!")
        println("$DASHEDDOUBLE${ANSI_RESET}")

    }

    // [VALIDATION] Minimal Nextflow version 
    //  We do this now and not on top because min_nf_version is a param itself that requires validation:
    if( !nextflow.version.matches(">=${params.min_nf_version}") ) {
        println "$ANSI_RED" + "$DASHEDDOUBLE"
        println "[VERSION ERROR] This workflow requires Nextflow version ${params.min_nf_version}"
        println "=> You are running version $nextflow.version."
        println "=> Use NXF_VER=${params.min_nf_version} nextflow run (...)"
        println "$DASHEDDOUBLE ${ANSI_RESET}"
        System.exit(1)
    }

    // Print params summary with adaptive spacing so columns are properly aligned regardless of param name
    // We use the same order as the schema params were defined in schema.nf for the summary report
    // [TODO] Maybe we print some boilderplate options on top of that in the future,
    //        such as container engine, GitHub URL etc.
    def max_char = params.keySet().collect { it.length() }.max()  
    println "$ANSI_GREEN" + "$DASHEDDOUBLE"
    println "[PARAMS SUMMARY]"
    println ""
    schema.each { name, entry -> 
        
        def use_length = max_char - name.length()
        def spacer = ' '.multiply(use_length)
        def m = entry['value']
        if(m!='') println "${name} ${spacer}:: ${m}" 

    }
    println "${DASHEDDOUBLE}"
    println "$ANSI_RESET"

}

/*
    This script is evaluated in main.nf so we run the function here and then return the params map
    so it is available in the global main.nf environment
*/
ValidateParams()
return(params)


