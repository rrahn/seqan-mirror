/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
  ============================================================================
  Copyright (C) 2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
  ============================================================================
  Author: Rene Maerker <rene.maerker@fu-berlin.de>
  12.01.2011
  ============================================================================
  Implements Listenerconcpet in seqan
  ==========================================================================*/

#ifndef LISTENER_BASIC_H_
#define LISTENER_BASIC_H_

namespace seqan {


	template <typename T>
	struct Listener
	{
		T container;

	};

	template <typename TListener, typename TValue>
	void
	notifyModified(TListener const & listener,
				   TValue * object)
	{

	}

}
#endif /* LISTENER_BASIC_H_ */
