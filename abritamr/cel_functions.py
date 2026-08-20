from cel import Context, evaluate


def evaluate_rule(rule:str, ctx: Context, name:str='row') -> bool:
    """
    Evaluate a CEL rule against the provided data.

    Parameters
    ----------
    rule : str
        The CEL rule to evaluate.
    data : dict
        Dictionary containing the data to be evaluated against the rule.
    name : str
        Name of the context variable.

    Returns
    -------
    bool
        Result of the rule evaluation (True or False).
    """

    try:
        result = evaluate(rule, ctx)
        return result
    except Exception as e:
        raise RuntimeError(f"Error evaluating rule: {e}")

    
def create_cel_context(data: dict, name:str='row') -> Context:

    """
    Create a CEL context with the provided data.

    Parameters
    ----------
    data : dict
        Dictionary containing the data to be added to the context.
    name : str
        Name of the context.

    Returns
    -------
    Context
        A CEL context populated with the provided data.
    """
    ctx = Context()
    ctx.add_variable(name, data)

    for cst in custom_rules():
        ctx.add_function(cst, custom_rules()[cst])

    return ctx


def custom_rules() -> dict:
    """
    Define custom CEL functions for use in rule evaluation.

    Returns
    -------
    dict
        Dictionary containing custom CEL functions.
    """
    return {
        "contains_any": contains_any,
    }

def contains_any(data:list, query:str) -> bool:

    """
    Custom CEL function to check if a class contains a query string.

    Parameters
    ----------
    data : list
        List of class strings.
    query : str
        Query string to search for.

    Returns
    -------
    bool
        True if the query string is found in any of the class strings, False otherwise.
    """
    for d in data:
        if query.lower() in d.lower():
            return True
    return False

