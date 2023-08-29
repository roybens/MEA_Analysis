
import maxlab.apicomm as comm


class WellPlate(object):
    """Interface to query various properties of the currently 
    configured wellplate.
    """
    
    @staticmethod
    def __send(command):
        with comm.api_context() as api:
            return api.send(command)
    
    def query_version(self):
        """Returns the integer for the wellplate version
        """
        return self.__send("wellplate_query_version")
    
    def query_rows(self):
        """Returns the number of well-rows in the wellplate
        """
        return self.__send("wellplate_query_rows")
    
    def query_columns(self):
        """Returns the number of well-rows in the wellplate
        """
        return self.__send("wellplate_query_columns")
