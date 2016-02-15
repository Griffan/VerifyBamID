/*
 *  Copyright (C) 2010  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __ERROR_HANDLER_H__
#define __ERROR_HANDLER_H__

#include <iostream>

/// Class that controls the handling of errors.
class ErrorHandler
{
public:

    /// This specifies how this class should respond to errors.
    enum HandlingType {EXCEPTION, ///< throw an exception for the error
                       ABORT,     ///< exit the program on the error
                       RETURN     ///< just return failure on the error
    };

    /// Constructor
    ErrorHandler();
   
    /// Destructor
    ~ErrorHandler();

    /// Handle an error based on the error handling type.
    static void handleError(const char* message, 
                            HandlingType handlingType = EXCEPTION);
      
private:
};


#endif
